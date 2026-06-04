"""生存时间 (HR)：临床预期 + 统计学设计 + 调度 + 结果。"""

import math as _math
import streamlit as st

import ui_common
from calc_survival import SurvivalCalculator


# ==============================================================================
# 模式选择器（表单外侧，瞬间响应）
# ==============================================================================
def _render_mode_selectors(
    group_type: str,
    test_type: str,
) -> dict:
    """渲染所有 st.radio / st.selectbox 模式选择器（在 st.form 外侧）。

    返回 dict 包含用户的选择，供后续决定表单内展示哪些数值输入框。
    """
    out: dict = {}

    if group_type == "单臂设计 (Single Arm)":
        hr_mode = st.radio(
            "HR 来源",
            ["由中位生存期推算（指数分布）", "直接输入 HR", "由生存率推算（如 1年 PFS 率）"],
            horizontal=True,
            key="surv_single_hr_mode",
        )
        if "生存率" in hr_mode:
            out["hr_input_mode"] = "survival_rate"
        elif "直接" in hr_mode:
            out["hr_input_mode"] = "direct"
        else:
            out["hr_input_mode"] = "median"

        ev_mode = st.radio(
            "事件概率",
            ["直接填写总体预期事件概率", "由入组期+随访推算（Lachin–Foulkes，指数）"],
            horizontal=True,
            key="surv_single_ev_mode",
        )
        out["event_mode"] = "lachin" if "入组" in ev_mode else "direct"

    else:
        # 独立两组 (Parallel)
        hr_mode = st.radio(
            "HR 来源",
            ["由中位生存期推算（指数分布）", "直接输入 HR", "由生存率推算（如 1年 PFS 率）"],
            horizontal=True,
            key="surv_two_hr_mode",
        )
        if "生存率" in hr_mode:
            out["hr_input_mode"] = "survival_rate"
        elif "直接" in hr_mode:
            out["hr_input_mode"] = "direct"
        else:
            out["hr_input_mode"] = "median"

        if test_type == "差异性检验":
            lr_lbl = st.selectbox(
                "Log-rank 样本量公式",
                ["Schoenfeld（常用）", "Freedman（小样本略保守）"],
                index=0,
                key="surv_two_lr_method",
            )
            out["lr_method"] = "freedman" if "Freedman" in lr_lbl else "schoenfeld"
        else:
            out["lr_method"] = "schoenfeld"

        ev_mode = st.radio(
            "总体事件概率",
            ["直接填写", "由入组期+随访推算（Lachin–Foulkes，指数；需两组中位生存期）"],
            horizontal=True,
            key="surv_two_ev_mode",
        )
        out["event_mode"] = "lachin" if "入组" in ev_mode else "direct"

    return out


# ==============================================================================
# 数值输入（表单内部）
# ==============================================================================
def _render_numeric_inputs(
    group_type: str,
    test_type: str,
    mode: dict,
) -> dict:
    """根据 mode 选择器状态，渲染对应的数值输入框（在 st.form 内部）。

    所有 st.number_input 均使用唯一 key，确保 Streamlit 正确跟踪。
    """
    out: dict = {}

    if group_type == "单臂设计 (Single Arm)":
        # ---- HR 模式 ----
        hr_input_mode = mode["hr_input_mode"]
        out["hr_input_mode"] = hr_input_mode

        if hr_input_mode == "survival_rate":
            t_m = st.number_input(
                "里程碑时间点 T（月）", value=12.0, step=0.5, min_value=0.1,
                key="surv_single_surv_t",
            )
            c1, c2 = st.columns(2)
            with c1:
                s_t_pct = st.number_input(
                    "试验组预期生存率 S_T(T) (%)",
                    value=50.0, step=1.0, min_value=0.1, max_value=99.9,
                    key="surv_single_surv_st",
                )
            with c2:
                s_ref_pct = st.number_input(
                    "参考生存率 S_0(T) (%)",
                    value=30.0, step=1.0, min_value=0.1, max_value=99.9,
                    key="surv_single_surv_sref",
                )
            out["survival_t"] = t_m
            out["s_t"] = s_t_pct / 100.0
            out["s_ref"] = s_ref_pct / 100.0
            out["median_t"] = None
            out["median_ref"] = None
            out["hr_direct"] = None
        elif hr_input_mode == "direct":
            out["hr_direct"] = st.number_input(
                "风险比 HR（试验/参考）",
                value=0.75, step=0.05, min_value=0.01, max_value=50.0,
                key="surv_single_hr_direct",
            )
            out["median_t"] = None
            out["median_ref"] = None
        else:
            c1, c2 = st.columns(2)
            with c1:
                out["median_t"] = st.number_input(
                    "试验组中位生存期（月）", value=12.0, step=0.5,
                    key="surv_single_med_t",
                )
            with c2:
                out["median_ref"] = st.number_input(
                    "参考中位生存期（月）", value=9.0, step=0.5,
                    key="surv_single_med_ref",
                )
            out["hr_direct"] = None

        # 界值
        if test_type == "非劣效性检验":
            out["margin_hr"] = st.number_input(
                "非劣界值 HR₀（须 >1）",
                value=1.3, step=0.05, min_value=1.001,
                key="surv_single_margin_hr",
            )
        else:
            out["margin_hr"] = None

        if test_type == "等效性检验":
            out["equiv_delta"] = st.number_input(
                "等效界值 Δ（HR 须落在 (1/Δ, Δ)）",
                value=1.25, step=0.05, min_value=1.001,
                key="surv_single_equiv_delta",
            )
        else:
            out["equiv_delta"] = None

        # 事件概率
        event_mode = mode["event_mode"]
        out["event_mode"] = event_mode
        if event_mode == "direct":
            out["event_rate_percent"] = st.number_input(
                "预计总体事件概率（%）",
                value=70.0, step=1.0, min_value=0.1, max_value=100.0,
                key="surv_single_ev_direct",
            )
            out["T_a"] = None
            out["T_f"] = None
        else:
            out["T_a"] = st.number_input(
                "入组期 T_a（月）", value=18.0, step=0.5, min_value=0.1,
                key="surv_single_ta",
            )
            out["T_f"] = st.number_input(
                "最短随访期 T_f（月）", value=12.0, step=0.5, min_value=0.0,
                key="surv_single_tf",
            )
            out["event_rate_percent"] = None
            if hr_input_mode == "direct":
                out["median_t_lf"] = st.number_input(
                    "试验组中位生存期（月，仅用于 Lachin–Foulkes 事件概率）",
                    value=12.0, step=0.5, min_value=0.1,
                    key="surv_single_med_t_lf",
                )
            else:
                out["median_t_lf"] = None

        return out

    # ==========================================================================
    # 独立两组 (Parallel)
    # ==========================================================================
    hr_input_mode = mode["hr_input_mode"]
    out["hr_input_mode"] = hr_input_mode

    if hr_input_mode == "survival_rate":
        t_m = st.number_input(
            "里程碑时间点 T（月）", value=12.0, step=0.5, min_value=0.1,
            key="surv_two_surv_t",
        )
        c1, c2 = st.columns(2)
        with c1:
            s_t_pct = st.number_input(
                "试验组预期生存率 S_T(T) (%)",
                value=50.0, step=1.0, min_value=0.1, max_value=99.9,
                key="surv_two_surv_st",
            )
        with c2:
            s_c_pct = st.number_input(
                "对照组预期生存率 S_C(T) (%)",
                value=30.0, step=1.0, min_value=0.1, max_value=99.9,
                key="surv_two_surv_sc",
            )
        out["survival_t"] = t_m
        out["s_t"] = s_t_pct / 100.0
        out["s_c"] = s_c_pct / 100.0
        out["median_t"] = None
        out["median_c"] = None
        out["hr_direct"] = None
    elif hr_input_mode == "direct":
        out["hr_direct"] = st.number_input(
            "风险比 HR（试验/对照）",
            value=0.75, step=0.05, min_value=0.01, max_value=50.0,
            key="surv_two_hr_direct",
        )
        out["median_t"] = None
        out["median_c"] = None
    else:
        c1, c2 = st.columns(2)
        with c1:
            out["median_t"] = st.number_input(
                "试验组中位生存期（月）", value=12.0, step=0.5,
                key="surv_two_med_t",
            )
        with c2:
            out["median_c"] = st.number_input(
                "对照组中位生存期（月）", value=9.0, step=0.5,
                key="surv_two_med_c",
            )
        out["hr_direct"] = None

    # 界值
    if test_type == "非劣效性检验":
        out["margin_hr"] = st.number_input(
            "非劣界值 HR₀（须 >1）",
            value=1.3, step=0.05, min_value=1.001,
            key="surv_two_margin_hr",
        )
    else:
        out["margin_hr"] = None

    if test_type == "等效性检验":
        out["equiv_delta"] = st.number_input(
            "等效界值 Δ（预期 HR 须在 (1/Δ, Δ) 内）",
            value=1.25, step=0.05, min_value=1.001,
            key="surv_two_equiv_delta",
        )
    else:
        out["equiv_delta"] = None

    # 事件概率
    event_mode = mode["event_mode"]
    out["lr_method"] = mode.get("lr_method", "schoenfeld")
    out["event_mode"] = event_mode
    if event_mode == "direct":
        out["event_rate_percent"] = st.number_input(
            "预计总体事件概率（%）",
            value=70.0, step=1.0, min_value=0.1, max_value=100.0,
            key="surv_two_ev_direct",
        )
        out["T_a"] = None
        out["T_f"] = None
        out["median_t_lf"] = None
        out["median_c_lf"] = None
    else:
        out["T_a"] = st.number_input(
            "入组期 T_a（月）", value=24.0, step=0.5, min_value=0.1,
            key="surv_two_ta",
        )
        out["T_f"] = st.number_input(
            "最短随访期 T_f（月）", value=12.0, step=0.5, min_value=0.0,
            key="surv_two_tf",
        )
        out["event_rate_percent"] = None
        if hr_input_mode == "direct":
            c1, c2 = st.columns(2)
            with c1:
                out["median_t_lf"] = st.number_input(
                    "试验组中位生存期（月，仅用于 Lachin）",
                    value=12.0, step=0.5, min_value=0.1,
                    key="surv_two_med_t_lf",
                )
            with c2:
                out["median_c_lf"] = st.number_input(
                    "对照组中位生存期（月，仅用于 Lachin）",
                    value=9.0, step=0.5, min_value=0.1,
                    key="surv_two_med_c_lf",
                )
            st.caption("直接输入 HR 时，两组中位仅用于推算各组随访期内事件概率，与 HR 无强制指数一致性。")
        else:
            out["median_t_lf"] = None
            out["median_c_lf"] = None
            st.caption("已用中位生存期推算 HR，并用于 Lachin–Foulkes 事件概率。")

    return out


# ==============================================================================
# 兼容旧接口
# ==============================================================================
def render_clinical_inputs(
    group_type: str,
    test_type: str = "优效性检验",
    *,
    legacy_simple: bool = False,
) -> dict:
    """兼容旧调用（如 app_legacy.py），一体式渲染所有控件。"""
    if legacy_simple:
        c_p1, c_p2 = st.columns(2)
        with c_p1:
            median_t = st.number_input("试验组中位生存期（月）", value=12.0, step=0.5)
        with c_p2:
            if group_type == "单臂设计 (Single Arm)":
                median_ref = st.number_input("参考中位生存期（月）", value=9.0, step=0.5)
                median_c = None
            else:
                median_c = st.number_input("对照组中位生存期（月）", value=9.0, step=0.5)
                median_ref = None
        event_rate_percent = st.number_input("预计事件发生率 (%)", value=70.0, step=1.0)
        return {
            "group_type": group_type, "test_type": test_type,
            "endpoint_mode": "hr", "hr_input_mode": "median",
            "hr_direct": None, "median_t": median_t,
            "median_ref": median_ref, "median_c": median_c,
            "event_rate_percent": event_rate_percent, "margin_percent": 0.0,
            "event_mode": "direct", "T_a": None, "T_f": None,
            "median_t_lf": None, "median_c_lf": None,
            "lr_method": "schoenfeld", "margin_hr": None, "equiv_delta": None,
        }

    # 非 legacy：一体式渲染所有控件（旧架构，仍可工作）
    mode = _render_mode_selectors(group_type, test_type)
    nums = _render_numeric_inputs(group_type, test_type, mode)
    return {**mode, **nums, "group_type": group_type, "test_type": test_type}


# ==============================================================================
# 主渲染入口
# ==============================================================================
def render(col_input, col_result, group_type: str, test_type: str) -> None:
    variable_type = "生存时间 (HR)"

    # 自定义提交按钮颜色
    st.markdown(
        """
        <style>
        div[data-testid="stFormSubmitButton"] button {
            background-color: #DE9F83 !important;
            border-color: #DE9F83 !important;
            color: white !important;
        }
        div[data-testid="stFormSubmitButton"] button:hover {
            background-color: #cb8d71 !important;
            border-color: #cb8d71 !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

    # ------------------------------------------------------------------
    # 左侧：模式选择器（表单外侧，瞬间响应）+ 数值输入（表单内部）
    # ------------------------------------------------------------------
    with col_input:
        st.subheader("临床预期", divider="gray")

        # 1. 模式选择器（表单外侧 — 切换即时响应）
        mode = _render_mode_selectors(group_type, test_type)

        # 2. 数值输入 + 统计学设计参数（表单内部）
        form_key = f"survival_form_{group_type}_{test_type}"
        with st.form(key=form_key):
            _s = _render_numeric_inputs(group_type, test_type, mode)

            # 统计学设计参数（也在表单内，与数值一并提交）
            stats = ui_common.render_statistics_design(variable_type, group_type, test_type)
            alpha = stats["alpha"]
            power_percent = stats["power_percent"]
            dropout_percent = stats["dropout_percent"]
            alloc_ratio_str = stats["alloc_ratio_str"]
            alpha_side_label = stats["alpha_side_label"]

            # ---- 提交按钮 ----
            st.markdown("<br>", unsafe_allow_html=True)
            submitted = st.form_submit_button(
                label="GO：开始计算样本量",
                type="primary",
                use_container_width=True,
            )

    # ------------------------------------------------------------------
    # 右侧：按钮未点击时显示提示；点击后执行计算并输出结果
    # ------------------------------------------------------------------
    if not submitted:
        with col_result:
            st.info(
                "👈 请在左侧填写临床预期与统计学设计参数，"
                "然后点击 **「GO：开始计算样本量」** 按钮。"
            )
        return

    # ------------------------------------------------------------------
    # 内部推导：由生存率推算中位生存期与 HR
    # ------------------------------------------------------------------
    hr_input_mode = _s.get("hr_input_mode", "median")
    if hr_input_mode == "survival_rate":
        t_m = _s["survival_t"]
        s_t = _s["s_t"]
        if group_type == "独立两组 (Parallel)":
            s_c = _s["s_c"]
            median_t = -t_m * _math.log(2) / _math.log(s_t)
            median_c = -t_m * _math.log(2) / _math.log(s_c)
            _s["median_t"] = median_t
            _s["median_c"] = median_c
        else:
            s_ref = _s["s_ref"]
            median_t = -t_m * _math.log(2) / _math.log(s_t)
            median_ref = -t_m * _math.log(2) / _math.log(s_ref)
            _s["median_t"] = median_t
            _s["median_ref"] = median_ref
        _s["hr_input_mode"] = "median"

    # ------------------------------------------------------------------
    # 计算
    # ------------------------------------------------------------------
    alloc_ratio_val = eval(alloc_ratio_str.replace(":", "/"))
    survival_calculator = SurvivalCalculator()
    error_message = None
    result = None
    try:
        if group_type == "独立两组 (Parallel)":
            _mt2 = _s.get("median_t") or _s.get("median_t_lf")
            _mc2 = _s.get("median_c") or _s.get("median_c_lf")
            result = survival_calculator.compute_two_arm(
                hr_input_mode=_s["hr_input_mode"],
                hr_direct=_s.get("hr_direct"),
                median_t=_mt2,
                median_c=_mc2,
                alpha=alpha,
                power=power_percent / 100,
                allocation_ratio=alloc_ratio_val,
                dropout_rate=dropout_percent / 100,
                test_type=test_type,
                margin_hr=_s.get("margin_hr"),
                equiv_delta=_s.get("equiv_delta"),
                event_mode=_s["event_mode"],
                event_rate=_s["event_rate_percent"] / 100.0
                if _s.get("event_rate_percent") is not None
                else None,
                T_a=_s.get("T_a"),
                T_f=_s.get("T_f"),
                lr_method=_s.get("lr_method", "schoenfeld"),
            )
        else:
            _mt = _s.get("median_t")
            if _mt is None and _s.get("median_t_lf") is not None:
                _mt = _s.get("median_t_lf")
            result = survival_calculator.compute_single_arm(
                hr_input_mode=_s["hr_input_mode"],
                hr_direct=_s.get("hr_direct"),
                median_t=_mt,
                median_ref=_s.get("median_ref"),
                alpha=alpha,
                power=power_percent / 100,
                dropout_rate=dropout_percent / 100,
                test_type=test_type,
                margin_hr=_s.get("margin_hr"),
                equiv_delta=_s.get("equiv_delta"),
                event_mode=_s["event_mode"],
                event_rate=_s["event_rate_percent"] / 100.0
                if _s.get("event_rate_percent") is not None
                else None,
                T_a=_s.get("T_a"),
                T_f=_s.get("T_f"),
            )
    except Exception as e:
        error_message = str(e)

    # ------------------------------------------------------------------
    # 右侧：结果输出
    # ------------------------------------------------------------------

    # 获取显示用的参数
    median_t_disp = _s.get("median_t") or _s.get("median_t_lf")
    median_ref_disp = _s.get("median_ref")
    median_c_disp = _s.get("median_c") or _s.get("median_c_lf")
    survival_t_disp = _s.get("survival_t")
    s_t_disp = _s.get("s_t")
    s_c_disp = _s.get("s_c")
    s_ref_disp = _s.get("s_ref")

    with col_result:
        if error_message:
            st.error(f"计算冲突: {error_message}")
            return
        if not result:
            st.warning("无计算结果。")
            return

        total_n = result["total_sample_size_with_dropout"]
        n_treat = result.get("n_treatment_with_dropout", 0)
        n_ctrl = result.get("n_control_with_dropout", 0)
        theo_total = result["total_sample_size"]

        if group_type == "单臂设计 (Single Arm)":
            st.markdown(
                f"""
<div class="result-kpi-card">
    <div style="font-size:16px;font-weight:500;text-align:center;opacity:0.95;letter-spacing:0.5px;">
        目标招募样本量
    </div>
    <div style="text-align:center;font-size:48px;font-weight:700;margin:24px 0;">
        {total_n} 例
    </div>
    <div class="total-box">
        理论所需：{theo_total} 例（已计入 {dropout_percent}% 脱落率）
    </div>
</div>
""",
                unsafe_allow_html=True,
            )

            st.markdown("<br>", unsafe_allow_html=True)
            st.markdown('<div class="section-title">样本量计算方法</div>', unsafe_allow_html=True)

            hr_show = result.get("hr", None)
            if survival_t_disp is not None and s_t_disp is not None and s_ref_disp is not None:
                cn_text = (
                    f"本研究为单臂生存终点{test_type[:2]}设计。"
                    f"根据临床预期，假定 {survival_t_disp:.1f} 月时试验组生存率为 {s_t_disp*100:.1f}%、"
                    f"参考生存率为 {s_ref_disp*100:.1f}%，"
                    f"推算风险比 HR≈{hr_show:.3f}，"
                    f"对应中位生存期约 {median_t_disp:.1f} 月与 {median_ref_disp:.1f} 月。"
                    f"采用{alpha_side_label}，显著性水平 α = {alpha}，检验效能 1-β = {power_percent:.0f}%。"
                    f"在上述参数条件下，依据 {result.get('method', '')} 计算，"
                    f"理论所需样本量为 {theo_total} 例。"
                    f"考虑约 {dropout_percent:.0f}% 脱落率，计划入组 {total_n} 例。"
                )
            elif median_t_disp is not None and median_ref_disp is not None:
                hr_clause = f"，HR≈{hr_show:.3f}" if hr_show is not None else ""
                cn_text = (
                    f"本研究为单臂生存终点{test_type[:2]}设计。"
                    f"试验组中位生存期 {median_t_disp:.1f} 月，"
                    f"参考中位生存期 {median_ref_disp:.1f} 月{hr_clause}。"
                    f"采用{alpha_side_label}，显著性水平 α = {alpha}，检验效能 1-β = {power_percent:.0f}%。"
                    f"在上述参数条件下，依据 {result.get('method', '')} 计算，"
                    f"理论所需样本量为 {theo_total} 例。"
                    f"考虑约 {dropout_percent:.0f}% 脱落率，计划入组 {total_n} 例。"
                )
            else:
                hr_clause = f"，HR={hr_show:.3f}" if hr_show is not None else ""
                cn_text = (
                    f"本研究为单臂生存终点{test_type[:2]}设计{hr_clause}。"
                    f"采用{alpha_side_label}，显著性水平 α = {alpha}，检验效能 1-β = {power_percent:.0f}%。"
                    f"在上述参数条件下，依据 {result.get('method', '')} 计算，"
                    f"理论所需样本量为 {theo_total} 例。"
                    f"考虑约 {dropout_percent:.0f}% 脱落率，计划入组 {total_n} 例。"
                )

            st.markdown(
                f"""
<div class="statement-box">
    <p>{cn_text}</p>
</div>
""",
                unsafe_allow_html=True,
            )
        else:
            # 独立两组：Morandi KPI 卡片
            st.markdown(
                f"""
<div class="result-kpi-card">
    <div style="font-size:16px;font-weight:500;text-align:center;opacity:0.95;letter-spacing:0.5px;">
        目标招募样本量
    </div>
    <div class="group-split">
        <div class="group-box">
            <div class="group-num">{n_treat}</div>
            <div class="group-name">试验组</div>
        </div>
        <div class="group-box">
            <div class="group-num">{n_ctrl}</div>
            <div class="group-name">对照组</div>
        </div>
    </div>
    <div class="total-box">
        总计招募：{total_n} 例（已计入 {dropout_percent}% 脱落率，分配比 {alloc_ratio_str}）
    </div>
</div>
""",
                unsafe_allow_html=True,
            )

            st.markdown("<br>", unsafe_allow_html=True)
            st.markdown('<div class="section-title">样本量计算方法</div>', unsafe_allow_html=True)

            hr_show = result.get("hr", None)
            if survival_t_disp is not None and s_t_disp is not None and s_c_disp is not None:
                cn_text = (
                    f"本研究为生存终点的独立两组{test_type[:2]}设计。"
                    f"根据临床预期，假定 {survival_t_disp:.1f} 月时试验组生存率为 {s_t_disp*100:.1f}%、"
                    f"对照组生存率为 {s_c_disp*100:.1f}%，"
                    f"推算风险比 HR≈{hr_show:.3f}，"
                    f"对应中位生存期约 {median_t_disp:.1f} 月与 {median_c_disp:.1f} 月。"
                    f"采用{alpha_side_label}，显著性水平 α = {alpha}，检验效能 1-β = {power_percent:.0f}%。"
                    f"在上述参数条件下，依据 {result.get('method', '')} 计算，"
                    f"理论总样本量为 {theo_total} 例。"
                    f"考虑约 {dropout_percent:.0f}% 脱落率，实际计划总入组 {total_n} 例"
                    f"（试验组 {n_treat} 例，对照组 {n_ctrl} 例）。"
                )
            elif median_t_disp is not None and median_c_disp is not None:
                hr_clause = f"，HR≈{hr_show:.3f}" if hr_show is not None else ""
                cn_text = (
                    f"本研究为生存终点的独立两组{test_type[:2]}设计。"
                    f"试验组中位生存期 {median_t_disp:.1f} 月，"
                    f"对照组中位生存期 {median_c_disp:.1f} 月{hr_clause}。"
                    f"采用{alpha_side_label}，显著性水平 α = {alpha}，检验效能 1-β = {power_percent:.0f}%。"
                    f"在上述参数条件下，依据 {result.get('method', '')} 计算，"
                    f"理论总样本量为 {theo_total} 例。"
                    f"考虑约 {dropout_percent:.0f}% 脱落率，实际计划总入组 {total_n} 例"
                    f"（试验组 {n_treat} 例，对照组 {n_ctrl} 例）。"
                )
            else:
                hr_clause = f"，HR={hr_show:.3f}" if hr_show is not None else ""
                cn_text = (
                    f"本研究为生存终点独立两组{test_type[:2]}设计{hr_clause}。"
                    f"采用{alpha_side_label}，显著性水平 α = {alpha}，检验效能 1-β = {power_percent:.0f}%。"
                    f"在上述参数条件下，依据 {result.get('method', '')} 计算，"
                    f"理论总样本量为 {theo_total} 例。"
                    f"考虑约 {dropout_percent:.0f}% 脱落率，实际计划总入组 {total_n} 例"
                    f"（试验组 {n_treat} 例，对照组 {n_ctrl} 例）。"
                )

            st.markdown(
                f"""
<div class="statement-box">
    <p>{cn_text}</p>
</div>
""",
                unsafe_allow_html=True,
            )
