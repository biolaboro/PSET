#!/usr/bin/env python3

import os
import sys
import shlex
import sqlite3
import time
from tempfile import NamedTemporaryFile
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from collections import Counter, namedtuple
from pathlib import Path
from subprocess import Popen, STDOUT


class Pool():

    Task = namedtuple("Task", ("row", "proc"))

    SQL = """
        BEGIN;
        CREATE TABLE IF NOT EXISTS task (
            "id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
            "user" TEXT NOT NULL,
            "args" TEXT NOT NULL,
            "cpu" INTEGER NOT NULL DEFAULT 1,
            "submitted" DATETIME DEFAULT (strftime('%Y-%m-%dT%H:%M:%S.%s', 'NOW')),
            "modified" DATETIME DEFAULT (strftime('%Y-%m-%dT%H:%M:%S.%s', 'NOW')),
            "status" TEXT NOT NULL CHECK (
                status == "QUEUED" OR
                status == "RUNNING" OR
                status == "SUCCESS" OR
                status == "FAILURE" OR
                status == "CANCELED" OR
                status == "DOCANCEL" OR
                status == "UNKNOWN"
            ),
            "log" TEXT DEFAULT ''
        );
        CREATE INDEX IF NOT EXISTS "task_status_idx" ON "task" ("status");
        CREATE TRIGGER IF NOT EXISTS task_modified_trg BEFORE UPDATE ON "task"
        BEGIN
            UPDATE "task"
            SET modified = strftime('%Y-%m-%dT%H:%M:%S:%s', 'now', 'localtime')
            WHERE id = old.id;
        END;
        COMMIT;
    """
    
    def __init__(self, db, cmd, wd, cpu=1, lag=1, lim=0, old=0):
        # args
        self.db, self.cmd, self.wd = db, cmd, wd
        # kwargs
        self.cpu, self.lag, self.lim, self.old = cpu, lag, lim, old
        # properties
        self.conn, self.pool, self.logs = None, set(), dict()

    def _connect(self):
        self.conn = sqlite3.connect(self.db)
        self.conn.row_factory = lambda curs, row: namedtuple("Row", (ele[0] for ele in curs.description))(*row)

    def _init_db(self):
        """initialize database connection and table"""
        # ensure table exists
        self.conn.executescript(self.SQL)
        # set all non-SUCCESS tasks to FAILURE
        self.conn.executescript("UPDATE task SET status = 'UNKNOWN' WHERE status != 'SUCCESS';")
        self.conn.commit()

    def _clean_task(self, ele, status):
        ## process log
        if value := self.logs.get(ele.row.id):
            with open(value.name) as log:
                sql = "UPDATE task SET status = :status, log = :log WHERE id == :id;"
                self.conn.execute(sql, dict(status=status, id=ele.row.id, log=log.read()))
                self.conn.commit()
            ## remove completed log from pool
            os.unlink(log.name)
            del self.logs[ele.row.id]
        ## wait (might not be necessary...)
        ele.proc.wait()
        ## remove completed task from pool
        if ele in self.pool:
            self.pool.remove(ele)

    def _process_old(self):
        if isinstance(days := self.old, int) and days > 0:
            self.conn.execute(f"DELETE FROM task WHERE modified < DATE('now', '-{days} days')")
            self.conn.commit()

    def _process_canceled(self):
        # if ids := set(row.id for row in self.conn.execute("SELECT id FROM task WHERE status == 'DOCANCEL';").fetchall()):
        ids = set(row.id for row in self.conn.execute("SELECT id FROM task WHERE status == 'DOCANCEL';").fetchall())
        # cancel active tasks
        if ids:
            for ele in list(self.pool):
                if ele.row.id in ids:
                    ele.proc.kill()
                    self._clean_task(ele, "CANCELED")
                    ids.remove(ele.row.id)
        # cancel tasks that never entered the pool
        self.conn.execute(f"UPDATE task SET status = 'CANCELED' WHERE status == 'DOCANCEL';")

    def _process_complete(self):
        """update status of completed tasks"""
        completed = {ele.row.id: ele for ele in self.completed()}
        sql = f"SELECT id FROM task WHERE id IN ({",".join("?" * len(completed))});"
        for row in self.conn.execute(sql, [ele.row.id for ele in completed.values()]):
            if ele := completed.get(row.id):
                status = "FAILURE" if ele.proc.returncode else "SUCCESS"
                self._clean_task(ele, status)
    
    def _process_running(self):
        """update status of running tasks"""
        for row in self.conn.execute("SELECT id, log FROM task WHERE status == 'RUNNING';").fetchall():
            with open(self.logs[row.id].name) as file:
                self.conn.execute("UPDATE task SET log = :log WHERE id == :id", dict(log=file.read(), id=row.id))
                self.conn.commit()
    
    def _process_queued(self):
        """update status of queued tasks"""
        # get the next tasks
        sql = "SELECT id, user, args, cpu, status FROM task WHERE status == 'QUEUED' ORDER BY submitted;"
        for row in self.conn.execute(sql):
            running_users = Counter(ele.row.user for ele in self.running())
            if row.cpu <= self.available() and (self.lim == 0 or running_users[row.user] <= self.lim):
                ## start process
                self.logs[row.id] = (log := NamedTemporaryFile(prefix=f"pool-{row.id}_", suffix=f".log", delete=False, delete_on_close=False))
                self.pool.add(self.Task(row, Popen((self.cmd, *shlex.split(row.args)), stdout=log, stderr=STDOUT, cwd=self.wd)))
                ## set RUNNING
                self.conn.execute("UPDATE task SET status = 'RUNNING' WHERE id == ?;", (row.id, ))
            else:   
                self.conn.execute("UPDATE task SET modified = strftime('%Y-%m-%dT%H:%M:%S:%s', 'now', 'localtime') WHERE id == ?;", (row.id, ))
            self.conn.commit()

    def available(self):
        """return the number of available CPUs"""
        return self.cpu - sum(ele.row.cpu for ele in self.running())

    def completed(self):
        """yield the next Task that is complete"""
        yield from (ele for ele in self.pool if ele.proc.poll() is not None)

    def running(self):
        """yield the next Task that is running"""
        yield from (ele for ele in self.pool if ele.proc.poll() is None)
    
    def status(self):
        return self.available(), self.cpu, list(self.conn.execute("SELECT * FROM TASK ORDER BY submitted;"))

    def status_print(self):
        """print status"""
        print(f"cpu: {self.available()} + {sum(1 for _ in self.running())} = {self.cpu}")
        for ele in self.pool:
            print(ele.row.id, ele.proc.returncode, self.logs[ele.row.id].name, os.path.exists(self.logs[ele.row.id].name))

    def run(self):
        """infinite run loop"""
        # connect/init
        self._connect()
        self._init_db()
        self.conn.close()
        # event loop
        while True:
            # self.status_print()
            self._connect()
            self._process_old()
            self._process_canceled()
            self._process_complete()
            self._process_running()
            self._process_queued()
            self.conn.close()
            time.sleep(self.lag)

def parse_args(args):
    parser = ArgumentParser(description="submission system", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("cmd", help="the command to execute")
    parser.add_argument("-wd", help="the working directory for execution", default=Path(__file__).parent.parent)
    parser.add_argument("-db", help="the path to the SQLite database", default=Path(__file__).parent / "pool.sdb")
    parser.add_argument("-cpu", help="the maximum number of CPUs", type=int, default=1)
    parser.add_argument("-lag", help="the event-loop lag in seconds", type=float, default=1)
    parser.add_argument("-lim", help="the maximum number of tasks per user (0 = unlimited)", type=int, default=1)
    parser.add_argument("-old", help="the max age of a task in the database (in days) (0 = unlimited)", type=int, default=30)
    return parser.parse_args(args)

def main(argv):
    args = parse_args(argv[1:])
    print(args)
    Pool(args.db, args.cmd, wd=args.wd, cpu=args.cpu, lag=args.lag, lim=args.lim, old=args.old).run()
    
if __name__ == "__main__":
    main(sys.argv)
