"""计算记录后台：密码进入后查看次数与明细。"""

from __future__ import annotations

import os
from datetime import datetime, timedelta

import pandas as pd
import streamlit as st

import usage_log


def _configured_password() -> str:
    env = os.environ.get("ADMIN_PASSWORD", "").strip()
    if env:
        return env
    try:
        return str(st.secrets["admin_password"]).strip()
    except Exception:
        return ""


def _hide_public_sidebar() -> None:
    st.markdown(
        """
        <style>
        [data-testid="stSidebar"],
        [data-testid="stSidebarCollapsedControl"] { display: none; }
        </style>
        """,
        unsafe_allow_html=True,
    )


def render() -> None:
    import ui_common

    ui_common.inject_morandi_ui()
    _hide_public_sidebar()
    password = _configured_password()
    if not password:
        st.markdown('<div class="main-title">计算记录后台</div>', unsafe_allow_html=True)
        st.warning(
            "后台密码还没设置。计算仍会写入记录，但这个页面先锁着。"
            "在服务器上设置环境变量 ADMIN_PASSWORD，"
            "或在 .streamlit/secrets.toml 里写 admin_password。"
        )
        return

    if not st.session_state.get("admin_authed"):
        st.markdown('<div class="main-title">计算记录后台</div>', unsafe_allow_html=True)
        _, mid, _ = st.columns([1, 1.1, 1])
        with mid:
            st.caption("内部页面。计算器里不会出现入口。")
            with st.form("admin_login"):
                entered = st.text_input("密码", type="password")
                submitted = st.form_submit_button("登录", type="primary", width="stretch")
            if submitted:
                if entered == password:
                    st.session_state["admin_authed"] = True
                    st.rerun()
                else:
                    st.error("密码不正确")
        return

    head_left, head_right = st.columns([5, 2])
    with head_left:
        st.markdown('<div class="main-title" style="text-align:left;margin-top:0;">计算记录后台</div>', unsafe_allow_html=True)
    with head_right:
        back, logout = st.columns(2)
        with back:
            st.link_button("返回计算器", "/", width="stretch")
        with logout:
            if st.button("退出", width="stretch"):
                st.session_state["admin_authed"] = False
                st.rerun()

    summary = usage_log.usage_summary()
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("累计计算", summary["total"])
    c2.metric("今日", summary["today"])
    c3.metric("近 7 日", summary["last_7_days"])
    c4.metric("近 30 日", summary["last_30_days"])

    if summary["total"] == 0:
        st.info("还没有计算记录。在计算器页面成功算出一次样本量后，这里会记上一笔。")
        return

    st.markdown('<div class="section-title">按时间段查询</div>', unsafe_allow_html=True)
    today = datetime.now().date()
    picked = st.date_input(
        "开始日期 — 结束日期",
        value=(today - timedelta(days=6), today),
        max_value=today,
    )
    if not (isinstance(picked, (tuple, list)) and len(picked) == 2):
        st.caption("请把开始日期和结束日期都选上。")
        return

    start_day, end_day = picked
    period = usage_log.query_range(
        datetime(start_day.year, start_day.month, start_day.day),
        datetime(end_day.year, end_day.month, end_day.day),
    )
    _render_period(period)
    _render_clear(summary["total"])


def _render_clear(total: int) -> None:
    st.markdown('<div class="section-title">清空记录</div>', unsafe_allow_html=True)
    if not st.session_state.get("confirm_clear"):
        st.caption("清空后累计、今日和各时间段的次数都会回到 0，不能恢复。")
        if st.button("清空全部记录"):
            st.session_state["confirm_clear"] = True
            st.rerun()
        return

    st.warning(f"将删除全部 {total} 条记录。")
    confirm, cancel = st.columns(2)
    with confirm:
        if st.button("确认清零", type="primary", width="stretch"):
            usage_log.clear_calculations()
            st.session_state["confirm_clear"] = False
            st.rerun()
    with cancel:
        if st.button("取消", width="stretch"):
            st.session_state["confirm_clear"] = False
            st.rerun()


def _busiest(rows: list[dict]) -> list[dict]:
    if not rows:
        return []
    top = max(row["n"] for row in rows)
    if top <= 0:
        return []
    return [row for row in rows if row["n"] == top]


def _hour_span(hour: str) -> str:
    start = int(hour)
    return f"{start:02d}:00–{(start + 1) % 24:02d}:00"


def _render_period(period: dict) -> None:
    peak_days = _busiest(period["by_day"])
    peak_hours = _busiest(period["by_hour"])
    m1, m2, m3 = st.columns(3)
    m1.metric("该时段次数", period["total"])
    if peak_days:
        day_label = "、".join(row["day"] for row in peak_days)
        m2.metric("最集中的日期", day_label, delta=f"{peak_days[0]['n']} 次", delta_color="off")
    else:
        m2.metric("最集中的日期", "—")
    if peak_hours:
        hour_label = "、".join(_hour_span(row["hour"]) for row in peak_hours)
        m3.metric("最集中的时段", hour_label, delta=f"{peak_hours[0]['n']} 次", delta_color="off")
    else:
        m3.metric("最集中的时段", "—")

    if period["total"] == 0:
        st.info("这个时间段没有计算。")
        return

    day_col, hour_col = st.columns(2)
    with day_col:
        st.markdown('<div class="section-title">按日</div>', unsafe_allow_html=True)
        day_df = pd.DataFrame(period["by_day"]).rename(columns={"day": "日期", "n": "次数"})
        st.bar_chart(day_df, x="日期", y="次数")
    with hour_col:
        st.markdown('<div class="section-title">一天中的集中时段</div>', unsafe_allow_html=True)
        hour_df = pd.DataFrame(
            [{"时段": _hour_span(row["hour"]), "次数": row["n"]} for row in period["by_hour"]]
        )
        st.bar_chart(hour_df, x="时段", y="次数")

    left, right = st.columns(2)
    with left:
        st.markdown('<div class="section-title">按终点变量</div>', unsafe_allow_html=True)
        st.dataframe(
            pd.DataFrame(period["by_variable"]).rename(columns={"name": "终点变量", "n": "次数"}),
            hide_index=True,
            width="stretch",
        )
        st.markdown('<div class="section-title">按检验目标</div>', unsafe_allow_html=True)
        st.dataframe(
            pd.DataFrame(period["by_test"]).rename(columns={"name": "检验目标", "n": "次数"}),
            hide_index=True,
            width="stretch",
        )
    with right:
        st.markdown('<div class="section-title">按分组设计</div>', unsafe_allow_html=True)
        st.dataframe(
            pd.DataFrame(period["by_group"]).rename(columns={"name": "分组设计", "n": "次数"}),
            hide_index=True,
            width="stretch",
        )

    st.markdown('<div class="section-title">该时段记录</div>', unsafe_allow_html=True)
    recent = pd.DataFrame(period["recent"]).rename(
        columns={
            "created_at": "时间",
            "variable_type": "终点变量",
            "group_type": "分组设计",
            "test_type": "检验目标",
        }
    )
    st.dataframe(recent, hide_index=True, width="stretch")
    if len(period["recent"]) == 50:
        st.caption("只列出该时段最近 50 条。上面的次数是全部记录。")
