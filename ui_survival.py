"""生存时间 (HR)：临床预期 + 统计学设计 + 调度 + 结果。"""

import streamlit as st

import ui_common
from calc_survival import SurvivalCalculator


def render_clinical_inputs(
    group_type: str,
    test_type: str = "优效性检验",
    *,
    legacy_simple: bool = False,
) -> dict:
    """单臂 / 双臂生存终点：HR 来源、界值、事件概率、Log-rank 方法。"""
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
            "group_type": group_type,
            "test_type": test_type,
            "endpoint_mode": "hr",
            "hr_input_mode": "median",
            "hr_direct": None,
            "median_t": median_t,
            "median_ref": median_ref,
            "median_c": median_c,
            "event_rate_percent": event_rate_percent,
            "margin_percent": 0.0,
            "event_mode": "direct",
            "T_a": None,
            "T_f": None,
            "median_t_lf": None,
            "median_c_lf": None,
            "lr_method": "schoenfeld",
            "margin_hr": None,
            "equiv_delta": None,
        }

    out: dict = {"group_type": group_type, "test_type": test_type}

    if group_type == "单臂设计 (Single Arm)":
        endpoint_mode = st.radio(
            "终点刻画方式",
            ["HR / 中位生存期", "里程碑生存率"],
            horizontal=True,
            help="里程碑：固定时间点 OS/PFS 率；HR：指数分布近似下相对历史对照。",
        )
        out["endpoint_mode"] = "milestone" if endpoint_mode.startswith("里程碑") else "hr"

        if out["endpoint_mode"] == "milestone":
            t_m = st.number_input("里程碑时间点（月）", value=12.0, step=0.5, min_value=0.1)
            c1, c2 = st.columns(2)
            with c1:
                s0_pct = st.number_input(
                    "参考生存率 S₀(t)（%）",
                    value=30.0,
                    step=1.0,
                    min_value=0.1,
                    max_value=99.9,
                    help="历史或外部对照在该时间点的生存率。",
                )
            with c2:
                s1_pct = st.number_input(
                    "试验组预期生存率 S₁(t)（%）",
                    value=45.0,
                    step=1.0,
                    min_value=0.1,
                    max_value=99.9,
                )
            out["milestone_t"] = t_m
            out["s0"] = s0_pct / 100.0
            out["s1"] = s1_pct / 100.0
            return out

        hr_mode = st.radio(
            "HR 来源",
            ["由中位生存期推算（指数分布）", "直接输入 HR"],
            horizontal=True,
        )
        out["hr_input_mode"] = "direct" if "直接" in hr_mode else "median"
        if out["hr_input_mode"] == "median":
            c1, c2 = st.columns(2)
            with c1:
                median_t = st.number_input("试验组中位生存期（月）", value=12.0, step=0.5)
            with c2:
                median_ref = st.number_input("参考中位生存期（月）", value=9.0, step=0.5)
            out["median_t"] = median_t
            out["median_ref"] = median_ref
            out["hr_direct"] = None
        else:
            out["hr_direct"] = st.number_input(
                "风险比 HR（试验/参考）",
                value=0.75,
                step=0.05,
                min_value=0.01,
                max_value=50.0,
                help="h_试验 / h_参考；试验获益时常 < 1。",
            )
            out["median_t"] = None
            out["median_ref"] = None

        if test_type == "非劣效性检验":
            out["margin_hr"] = st.number_input(
                "非劣界值 HR₀（须 >1）",
                value=1.3,
                step=0.05,
                min_value=1.001,
                help="需证明试验 HR < HR₀。",
            )
        else:
            out["margin_hr"] = None

        if test_type == "等效性检验":
            out["equiv_delta"] = st.number_input(
                "等效界值 Δ（HR 须落在 (1/Δ, Δ)）",
                value=1.25,
                step=0.05,
                min_value=1.001,
            )
        else:
            out["equiv_delta"] = None

        ev_mode = st.radio(
            "事件概率",
            ["直接填写总体预期事件概率", "由入组期+随访推算（Lachin–Foulkes，指数）"],
            horizontal=True,
        )
        out["event_mode"] = "lachin" if "入组" in ev_mode else "direct"
        if out["event_mode"] == "direct":
            out["event_rate_percent"] = st.number_input(
                "预计总体事件概率（%）",
                value=70.0,
                step=1.0,
                min_value=0.1,
                max_value=100.0,
            )
            out["T_a"] = None
            out["T_f"] = None
        else:
            out["T_a"] = st.number_input("入组期 T_a（月）", value=18.0, step=0.5, min_value=0.1)
            out["T_f"] = st.number_input("最短随访期 T_f（月）", value=12.0, step=0.5, min_value=0.0)
            out["event_rate_percent"] = None
            if out["hr_input_mode"] == "direct":
                out["median_t_lf"] = st.number_input(
                    "试验组中位生存期（月，仅用于 Lachin–Foulkes 事件概率）",
                    value=12.0,
                    step=0.5,
                    min_value=0.1,
                    help="直接输入 HR 时需单独给出中位生存期以推算随访期内事件概率。",
                )
            else:
                out["median_t_lf"] = None

        return out

    # 双臂
    hr_mode = st.radio(
        "HR 来源",
        ["由中位生存期推算（指数分布）", "直接输入 HR"],
        horizontal=True,
    )
    out["hr_input_mode"] = "direct" if "直接" in hr_mode else "median"
    if out["hr_input_mode"] == "median":
        c1, c2 = st.columns(2)
        with c1:
            median_t = st.number_input("试验组中位生存期（月）", value=12.0, step=0.5)
        with c2:
            median_c = st.number_input("对照组中位生存期（月）", value=9.0, step=0.5)
        out["median_t"] = median_t
        out["median_c"] = median_c
        out["hr_direct"] = None
    else:
        out["hr_direct"] = st.number_input(
            "风险比 HR（试验/对照）",
            value=0.75,
            step=0.05,
            min_value=0.01,
            max_value=50.0,
        )
        out["median_t"] = None
        out["median_c"] = None

    if test_type == "非劣效性检验":
        out["margin_hr"] = st.number_input(
            "非劣界值 HR₀（须 >1）",
            value=1.3,
            step=0.05,
            min_value=1.001,
        )
    else:
        out["margin_hr"] = None

    if test_type == "等效性检验":
        out["equiv_delta"] = st.number_input(
            "等效界值 Δ（预期 HR 须在 (1/Δ, Δ) 内）",
            value=1.25,
            step=0.05,
            min_value=1.001,
        )
    else:
        out["equiv_delta"] = None

    if test_type == "差异性检验":
        lr_lbl = st.selectbox(
            "Log-rank 样本量公式",
            ["Schoenfeld（常用）", "Freedman（小样本略保守）"],
            index=0,
        )
        out["lr_method"] = "freedman" if "Freedman" in lr_lbl else "schoenfeld"
    else:
        out["lr_method"] = "schoenfeld"

    ev_mode = st.radio(
        "总体事件概率",
        ["直接填写", "由入组期+随访推算（Lachin–Foulkes，指数；需两组中位生存期）"],
        horizontal=True,
    )
    out["event_mode"] = "lachin" if "入组" in ev_mode else "direct"
    if out["event_mode"] == "direct":
        out["event_rate_percent"] = st.number_input(
            "预计总体事件概率（%）",
            value=70.0,
            step=1.0,
            min_value=0.1,
            max_value=100.0,
        )
        out["T_a"] = None
        out["T_f"] = None
        out["median_t_lf"] = None
        out["median_c_lf"] = None
    else:
        out["T_a"] = st.number_input("入组期 T_a（月）", value=24.0, step=0.5, min_value=0.1)
        out["T_f"] = st.number_input("最短随访期 T_f（月）", value=12.0, step=0.5, min_value=0.0)
        out["event_rate_percent"] = None
        if out["hr_input_mode"] == "direct":
            c1, c2 = st.columns(2)
            with c1:
                out["median_t_lf"] = st.number_input(
                    "试验组中位生存期（月，仅用于 Lachin）",
                    value=12.0,
                    step=0.5,
                    min_value=0.1,
                )
            with c2:
                out["median_c_lf"] = st.number_input(
                    "对照组中位生存期（月，仅用于 Lachin）",
                    value=9.0,
                    step=0.5,
                    min_value=0.1,
                )
            st.caption("直接输入 HR 时，两组中位仅用于推算各组随访期内事件概率，与 HR 无强制指数一致性。")
        else:
            out["median_t_lf"] = None
            out["median_c_lf"] = None
            st.caption("已用中位生存期推算 HR，并用于 Lachin–Foulkes 事件概率。")

    return out


def render(col_input, col_result, group_type: str, test_type: str) -> None:
    variable_type = "生存时间 (HR)"
    with col_input:
        st.markdown('<div class="section-title">临床预期</div>', unsafe_allow_html=True)
        _s = render_clinical_inputs(group_type, test_type)
        stats = ui_common.render_statistics_design(variable_type, group_type, test_type)
        alpha = stats["alpha"]
        power_percent = stats["power_percent"]
        dropout_percent = stats["dropout_percent"]
        alloc_ratio_str = stats["alloc_ratio_str"]
        alpha_side_label = stats["alpha_side_label"]

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
        elif _s.get("endpoint_mode") == "milestone":
            result = survival_calculator.compute_single_arm_milestone(
                t_months=_s["milestone_t"],
                s0=_s["s0"],
                s1=_s["s1"],
                alpha=alpha,
                power=power_percent / 100,
                dropout_rate=dropout_percent / 100,
                test_type=test_type,
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

    median_t_disp = _s.get("median_t") or _s.get("median_t_lf")
    median_ref_disp = _s.get("median_ref")
    median_c_disp = _s.get("median_c") or _s.get("median_c_lf")
    event_rate_pct_disp = _s.get("event_rate_percent")
    endpoint_mode = _s.get("endpoint_mode", "hr")

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
        theo_treat = int(result.get("n_treatment", theo_total / (1 + 1 / alloc_ratio_val)))
        theo_ctrl = int(result.get("n_control", theo_total - theo_treat))

        ev_p = result.get("event_prob_overall")
        ev_note = ""
        if ev_p is not None:
            ev_note = f"推算/假定总体事件概率≈{ev_p*100:.1f}%。"

        if group_type == "单臂设计 (Single Arm)":
            st.markdown(
                f"""
                <div class="result-kpi-card">
                    <div style="font-size: 16px; font-weight: 500; text-align: center; opacity: 0.95;">目标招募样本量</div>
                    <div style="text-align: center; margin: 30px 0;">
                        <span style="font-size: 64px; font-weight: 700;">{total_n}</span>
                    </div>
                    <div class="total-box">理论所需：{theo_total} 例 <span style="font-weight: 400; font-size: 15px; opacity: 0.9;">(已计入 {dropout_percent}% 脱落率)</span></div>
                </div>
                """,
                unsafe_allow_html=True,
            )
            st.markdown('<div class="section-title" style="margin-top:0;">样本量计算方法</div>', unsafe_allow_html=True)
            _meth = result.get("method", "")
            if endpoint_mode == "milestone":
                s0 = result.get("s0", _s.get("s0"))
                s1 = result.get("s1", _s.get("s1"))
                tm = result.get("milestone_t", _s.get("milestone_t"))
                cn_text_p1 = (
                    f"本研究为单臂生存终点**里程碑**分析（{tm:.1f} 月）。参考生存率 S₀≈{s0*100:.1f}%，"
                    f"预期 S₁≈{s1*100:.1f}%。采用{alpha_side_label}检验，"
                    f"显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                )
                cn_text_p2 = (
                    f"样本量采用【{_meth}】。理论所需样本量为 {theo_total} 例；"
                    f"考虑约 {dropout_percent:.1f}% 脱落率，计划入组 {total_n} 例。"
                )
            else:
                hr_show = result.get("hr", None)
                event_show = result.get("events", None)
                hr_clause = f"，对应预期 HR≈{hr_show:.3f}" if hr_show is not None else ""
                evt_clause = f"预计所需事件数约 {event_show} 例。" if event_show is not None else ""
                er_clause = (
                    f"预计总体事件概率为 {event_rate_pct_disp:.1f}% 时，"
                    if event_rate_pct_disp is not None
                    else f"{ev_note} "
                )
                if median_t_disp is not None and median_ref_disp is not None:
                    cn_text_p1 = (
                        f"本研究为单臂生存终点（相对历史对照）{test_type[:2]}设计。假定试验组中位生存期 {median_t_disp:.1f} 月，"
                        f"参考中位生存期 {median_ref_disp:.1f} 月{hr_clause}。采用{alpha_side_label}检验，"
                        f"显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                    )
                else:
                    cn_text_p1 = (
                        f"本研究为单臂生存终点{test_type[:2]}设计（HR 直接输入或中位生存期部分缺失）。"
                        f"{hr_clause} 采用{alpha_side_label}检验，显著性水平 α = {alpha}，"
                        f"检验效能（1−β）= {power_percent:.1f}%。"
                    )
                cn_text_p2 = (
                    f"{evt_clause}{er_clause}理论所需样本量为 {theo_total} 例。"
                    f"考虑约 {dropout_percent:.1f}% 脱落率，计划入组 {total_n} 例。"
                )
            st.markdown(
                f"<div class='statement-box'><p>{cn_text_p1}</p><p>{cn_text_p2}</p></div>",
                unsafe_allow_html=True,
            )
        else:
            st.markdown(
                f"""
                <div class="result-kpi-card">
                    <div style="font-size: 16px; font-weight: 500; text-align: center; opacity: 0.95; letter-spacing: 0.5px;">目标招募样本量</div>
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
                        总计招募：{total_n} 例 <span style="font-weight: 400; font-size: 15px; opacity: 0.9;">(已计入 {dropout_percent}% 脱落率)</span>
                    </div>
                </div>
                """,
                unsafe_allow_html=True,
            )
            st.markdown('<div class="section-title" style="margin-top:0;">样本量计算方法</div>', unsafe_allow_html=True)
            alloc_txt = f"即每组 {theo_treat} 例" if alloc_ratio_val == 1.0 else f"即试验组 {theo_treat} 例，对照组 {theo_ctrl} 例"
            hr_show = result.get("hr", None)
            event_show = result.get("events", None)
            hr_clause = f"，对应预期 HR≈{hr_show:.3f}" if hr_show is not None else ""
            evt_clause = f"所需事件数约 {event_show} 例，" if event_show is not None else ""
            er_txt = (
                f"预计总体事件概率 {event_rate_pct_disp:.1f}% 时，"
                if event_rate_pct_disp is not None
                else ev_note
            )
            if median_t_disp is not None and median_c_disp is not None:
                cn_text_p1 = (
                    f"本研究为生存终点的独立两组{test_type[:2]}设计。假定试验组中位生存期 {median_t_disp:.1f} 月，"
                    f"对照组中位生存期 {median_c_disp:.1f} 月{hr_clause}。采用{alpha_side_label}检验，"
                    f"显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                )
            else:
                cn_text_p1 = (
                    f"本研究为生存终点独立两组{test_type[:2]}设计（HR 直接输入）。"
                    f"{hr_clause} 采用{alpha_side_label}检验，显著性水平 α = {alpha}，"
                    f"检验效能（1−β）= {power_percent:.1f}%。"
                )
            _meth = result.get("method", "")
            cn_text_p2 = (
                f"样本量采用【{_meth}】。{er_txt}{evt_clause}理论总样本量为 {theo_total} 例，{alloc_txt}。"
                f"考虑约 {dropout_percent:.1f}% 脱落率，实际计划总入组 {total_n} 例，"
                f"其中试验组 {n_treat} 例，对照组 {n_ctrl} 例。"
            )
            st.markdown(
                f"<div class='statement-box'><p>{cn_text_p1}</p><p>{cn_text_p2}</p></div>",
                unsafe_allow_html=True,
            )
