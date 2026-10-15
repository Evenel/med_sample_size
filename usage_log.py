"""计算次数记录：成功算出样本量时写入 SQLite，供后台汇总。"""

from __future__ import annotations

import os
import sqlite3
import sys
import threading
from datetime import datetime, timedelta, timezone
from pathlib import Path

_LOCK = threading.Lock()
_TZ = timezone(timedelta(hours=8))


def db_path() -> Path:
    override = os.environ.get("CALC_LOG_PATH", "").strip()
    if override:
        return Path(override)
    return Path(__file__).resolve().parent / "data" / "calc_log.sqlite"


def _connect() -> sqlite3.Connection:
    path = db_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(path, timeout=10)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS calc_log (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            created_at TEXT NOT NULL,
            variable_type TEXT NOT NULL,
            group_type TEXT NOT NULL,
            test_type TEXT NOT NULL
        )
        """
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_calc_log_created_at ON calc_log (created_at)"
    )
    return conn


def _now() -> datetime:
    return datetime.now(_TZ)


def clear_calculations() -> int:
    """删除全部计算记录，次数回到 0。返回删除条数。"""
    with _LOCK:
        conn = _connect()
        try:
            deleted = int(conn.execute("SELECT COUNT(*) AS n FROM calc_log").fetchone()["n"])
            conn.execute("DELETE FROM calc_log")
            try:
                conn.execute("DELETE FROM sqlite_sequence WHERE name = 'calc_log'")
            except sqlite3.OperationalError:
                pass
            conn.commit()
            return deleted
        finally:
            conn.close()


def record_calculation(variable_type: str, group_type: str, test_type: str) -> None:
    """记一次成功计算。写入失败时不影响页面出结果。"""
    created_at = _now().isoformat(timespec="seconds")
    try:
        with _LOCK:
            conn = _connect()
            try:
                conn.execute(
                    """
                    INSERT INTO calc_log (created_at, variable_type, group_type, test_type)
                    VALUES (?, ?, ?, ?)
                    """,
                    (created_at, variable_type, group_type, test_type),
                )
                conn.commit()
            finally:
                conn.close()
    except Exception as exc:
        print(f"calc log write failed: {exc}", file=sys.stderr)


def _count_since(conn: sqlite3.Connection, start_iso: str) -> int:
    row = conn.execute(
        "SELECT COUNT(*) AS n FROM calc_log WHERE created_at >= ?",
        (start_iso,),
    ).fetchone()
    return int(row["n"])


def _group_counts(
    conn: sqlite3.Connection,
    column: str,
    start_iso: str | None = None,
    end_iso: str | None = None,
) -> list[dict]:
    if column not in {"variable_type", "group_type", "test_type"}:
        raise ValueError(column)
    sql = f"""
        SELECT {column} AS name, COUNT(*) AS n
        FROM calc_log
    """
    params: list[str] = []
    if start_iso is not None and end_iso is not None:
        sql += " WHERE created_at >= ? AND created_at < ?"
        params.extend([start_iso, end_iso])
    sql += f" GROUP BY {column} ORDER BY n DESC, name ASC"
    rows = conn.execute(sql, params).fetchall()
    return [{"name": row["name"], "n": int(row["n"])} for row in rows]


def usage_summary(recent_limit: int = 50) -> dict:
    now = _now()
    today_start = now.replace(hour=0, minute=0, second=0, microsecond=0)
    week_start = today_start - timedelta(days=6)
    month_start = today_start - timedelta(days=29)
    with _LOCK:
        conn = _connect()
        try:
            total = int(conn.execute("SELECT COUNT(*) AS n FROM calc_log").fetchone()["n"])
            by_day_rows = conn.execute(
                """
                SELECT substr(created_at, 1, 10) AS day, COUNT(*) AS n
                FROM calc_log
                WHERE created_at >= ?
                GROUP BY day
                ORDER BY day ASC
                """,
                (month_start.isoformat(timespec="seconds"),),
            ).fetchall()
            recent_rows = conn.execute(
                """
                SELECT created_at, variable_type, group_type, test_type
                FROM calc_log
                ORDER BY id DESC
                LIMIT ?
                """,
                (recent_limit,),
            ).fetchall()
            return {
                "total": total,
                "today": _count_since(conn, today_start.isoformat(timespec="seconds")),
                "last_7_days": _count_since(conn, week_start.isoformat(timespec="seconds")),
                "last_30_days": _count_since(conn, month_start.isoformat(timespec="seconds")),
                "by_variable": _group_counts(conn, "variable_type"),
                "by_group": _group_counts(conn, "group_type"),
                "by_test": _group_counts(conn, "test_type"),
                "by_day": [{"day": row["day"], "n": int(row["n"])} for row in by_day_rows],
                "recent": [dict(row) for row in recent_rows],
            }
        finally:
            conn.close()


def _midnight(value: datetime) -> datetime:
    """日历日按东八区的 00:00。无时区的日期视为东八区当天，不按机器时区换算。"""
    day = value.date() if value.tzinfo is None else value.astimezone(_TZ).date()
    return datetime(day.year, day.month, day.day, tzinfo=_TZ)


def _bounds(start: datetime, end: datetime) -> tuple[str, str]:
    """开始日 00:00 含，结束日次日 00:00 不含。"""
    start_dt = _midnight(start)
    end_next = _midnight(end) + timedelta(days=1)
    return start_dt.isoformat(timespec="seconds"), end_next.isoformat(timespec="seconds")


def query_range(start: datetime, end: datetime, recent_limit: int = 50) -> dict:
    """某个时间段内的次数、按日、按小时，以及该时段内的分类汇总。"""
    if end < start:
        start, end = end, start
    start_iso, end_iso = _bounds(start, end)
    with _LOCK:
        conn = _connect()
        try:
            total = int(
                conn.execute(
                    "SELECT COUNT(*) AS n FROM calc_log WHERE created_at >= ? AND created_at < ?",
                    (start_iso, end_iso),
                ).fetchone()["n"]
            )
            by_day_rows = conn.execute(
                """
                SELECT substr(created_at, 1, 10) AS day, COUNT(*) AS n
                FROM calc_log
                WHERE created_at >= ? AND created_at < ?
                GROUP BY day
                ORDER BY day ASC
                """,
                (start_iso, end_iso),
            ).fetchall()
            by_hour_rows = conn.execute(
                """
                SELECT substr(created_at, 12, 2) AS hour, COUNT(*) AS n
                FROM calc_log
                WHERE created_at >= ? AND created_at < ?
                GROUP BY hour
                ORDER BY hour ASC
                """,
                (start_iso, end_iso),
            ).fetchall()
            recent_rows = conn.execute(
                """
                SELECT created_at, variable_type, group_type, test_type
                FROM calc_log
                WHERE created_at >= ? AND created_at < ?
                ORDER BY id DESC
                LIMIT ?
                """,
                (start_iso, end_iso, recent_limit),
            ).fetchall()
            day_counts = {row["day"]: int(row["n"]) for row in by_day_rows}
            hour_counts = {row["hour"]: int(row["n"]) for row in by_hour_rows}
            start_day = datetime.fromisoformat(start_iso).date()
            end_day = datetime.fromisoformat(end_iso).date() - timedelta(days=1)
            by_day = []
            day = start_day
            while day <= end_day:
                key = day.isoformat()
                by_day.append({"day": key, "n": day_counts.get(key, 0)})
                day += timedelta(days=1)
            by_hour = [
                {"hour": f"{h:02d}", "n": hour_counts.get(f"{h:02d}", 0)}
                for h in range(24)
            ]
            return {
                "total": total,
                "by_day": by_day,
                "by_hour": by_hour,
                "by_variable": _group_counts(conn, "variable_type", start_iso, end_iso),
                "by_group": _group_counts(conn, "group_type", start_iso, end_iso),
                "by_test": _group_counts(conn, "test_type", start_iso, end_iso),
                "recent": [dict(row) for row in recent_rows],
            }
        finally:
            conn.close()
