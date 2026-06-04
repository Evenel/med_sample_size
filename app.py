"""临床试验样本量计算 — 入口：侧边栏 + 路由到各变量类型模块。"""

import streamlit as st

import ui_common
import ui_continuous
import ui_proportion
import ui_survival

st.set_page_config(
    page_title="临床试验样本量计算",
    layout="wide",
    initial_sidebar_state="expanded",
)

_SUPPORTED = {
    ("分类变量 (率)", "独立两组 (Parallel)"),
    ("分类变量 (率)", "单臂设计 (Single Arm)"),
    ("分类变量 (率)", "Simon 二阶段"),
    ("分类变量 (率)", "三组设计 (Three Arm)"),
    ("连续变量 (均值)", "独立两组 (Parallel)"),
    ("连续变量 (均值)", "单臂设计 (Single Arm)"),
    ("连续变量 (均值)", "配对两组 (Paired)"),
    ("连续变量 (均值)", "三组设计 (Three Arm)"),
    ("分类变量 (率)", "配对两组 (Paired)"),
    ("生存时间 (HR)", "独立两组 (Parallel)"),
    ("生存时间 (HR)", "单臂设计 (Single Arm)"),
}


def main_layout(variable_type, group_type, test_type, simon_design_type):
    ui_common.inject_morandi_ui()
    st.markdown(
        '<div class="main-title">临床试验样本量计算</div>',
        unsafe_allow_html=True,
    )
    if (variable_type, group_type) not in _SUPPORTED:
        st.warning(
            f"当前组合【{variable_type} + {group_type}】已纳入菜单，但算法模块仍在开发中。"
        )
        return

    col_input, col_result = st.columns([5, 5], gap="large")
    if variable_type == "分类变量 (率)":
        ui_proportion.render(col_input, col_result, group_type, test_type, simon_design_type)
    elif variable_type == "连续变量 (均值)":
        ui_continuous.render(col_input, col_result, group_type, test_type)
    else:
        ui_survival.render(col_input, col_result, group_type, test_type)


if __name__ == "__main__":
    v, g, t, sd = ui_common.sidebar_layout()
    main_layout(v, g, t, sd)
