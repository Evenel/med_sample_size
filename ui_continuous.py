"""连续变量 (均值)：临床预期 + 统计学设计 + 调度 + 结果。"""

import math

import streamlit as st

import ui_common
from calc_continuous import ContinuousCalculator


def render_clinical_inputs(group_type: str, test_type: str) -> dict:
    """双臂 / 单臂 / 配对均值的临床预期。"""
    if group_type == "独立两组 (Parallel)":
        if test_type == "精度分析":
            c_p1, c_p2 = st.columns(2)
            with c_p1:
                mean_t = st.number_input("试验组预期均值（描述用）", value=10.0, step=0.1)
            with c_p2:
                mean_c = st.number_input("对照组预期均值（描述用）", value=8.0, step=0.1)
            c_sd1, c_sd2 = st.columns(2)
            with c_sd1:
                sd_t = st.number_input("试验组标准差 SD", value=3.0, step=0.1)
            with c_sd2:
                sd_c = st.number_input("对照组标准差 SD", value=3.0, step=0.1)
            half_width_abs = st.number_input(
                "均值差 (μT−μC) 的 CI 半宽 w（与均值同单位）",
                value=0.5,
                step=0.05,
                min_value=0.001,
                help="双侧 (1−α) 置信区间对均值差的半宽须 ≤ w；按 Welch SE + t 分位数迭代。",
            )
            return {
                "mean_t": mean_t,
                "mean_c": mean_c,
                "sd_t": sd_t,
                "sd_c": sd_c,
                "margin_percent": 0.0,
                "half_width_abs": half_width_abs,
                "precision_mode": True,
                "two_arm_precision": True,
            }
        c_p1, c_p2 = st.columns(2)
        with c_p1:
            mean_t = st.number_input("试验组预期均值", value=10.0, step=0.1)
        with c_p2:
            mean_c = st.number_input("对照组预期均值", value=8.0, step=0.1)
        c_sd1, c_sd2 = st.columns(2)
        with c_sd1:
            sd_t = st.number_input("试验组标准差 SD", value=3.0, step=0.1)
        with c_sd2:
            sd_c = st.number_input("对照组标准差 SD", value=3.0, step=0.1)
        margin_percent = 0.0 if test_type == "差异性检验" else st.number_input("界值 Δ（绝对均值差）", value=1.0, step=0.1)
        two_arm_method = "pooled_t"
        if test_type == "差异性检验":
            _tw = st.selectbox(
                "计算方法（差异性）",
                [
                    "两独立样本 t（合并方差）",
                    "Welch t（方差不齐）",
                    "Wilcoxon 秩和 / Mann-Whitney（ARE 校正）",
                ],
                index=0,
                help="优效/非劣/等效用两样本 t + 界值；等效用 TOST。",
            )
            _twmap = {
                "两独立样本 t（合并方差）": "pooled_t",
                "Welch t（方差不齐）": "welch",
                "Wilcoxon 秩和 / Mann-Whitney（ARE 校正）": "wilcoxon_mw",
            }
            two_arm_method = _twmap[_tw]
        return {
            "mean_t": mean_t,
            "mean_c": mean_c,
            "sd_t": sd_t,
            "sd_c": sd_c,
            "margin_percent": margin_percent,
            "two_arm_method": two_arm_method,
            "precision_mode": False,
            "two_arm_precision": False,
        }
    if group_type == "单臂设计 (Single Arm)":
        if test_type == "精度分析":
            mean_t = st.number_input("预期均值 μ（用于方案描述，可选）", value=10.0, step=0.1)
            mean_ref = mean_t
            sd_single = st.number_input("预期总体标准差 SD", value=3.0, step=0.1)
            half_width_abs = st.number_input(
                "均值 CI 半宽 w（与均值同单位）",
                value=0.5,
                step=0.05,
                min_value=0.001,
                help="双侧 (1−α) 置信区间半宽须 ≤ w；样本量按 t 分布迭代。",
            )
            return {
                "mean_t": mean_t,
                "mean_ref": mean_ref,
                "sd_single": sd_single,
                "margin_percent": 0.0,
                "half_width_abs": half_width_abs,
                "precision_mode": True,
            }
        c_p1, c_p2 = st.columns(2)
        with c_p1:
            mean_t = st.number_input("试验组预期均值", value=10.0, step=0.1)
        with c_p2:
            mean_ref = st.number_input("参考均值（历史/目标）", value=8.0, step=0.1)
        sd_single = st.number_input("总体标准差 SD", value=3.0, step=0.1)
        margin_percent = 0.0 if test_type == "差异性检验" else st.number_input("界值 Δ（绝对均值差）", value=1.0, step=0.1)
        single_arm_method = "t"
        if test_type == "差异性检验":
            _lbl = st.selectbox(
                "计算方法（差异性）",
                ["单样本 t 检验", "Wilcoxon 符号秩（ARE 校正）"],
                index=0,
                help="优效/非劣/等效固定为单样本 t；差异性可在 t 与 Wilcoxon 间选择。",
            )
            single_arm_method = "wilcoxon" if _lbl.startswith("Wilcoxon") else "t"
        return {
            "mean_t": mean_t,
            "mean_ref": mean_ref,
            "sd_single": sd_single,
            "margin_percent": margin_percent,
            "single_arm_method": single_arm_method,
            "precision_mode": False,
        }
    if group_type == "配对两组 (Paired)":
        sd_t = None
        sd_c = None
        rho = None
        sd_diff = st.number_input(
            "配对差值标准差 SD_diff",
            value=3.0,
            step=0.1,
            min_value=1e-6,
        )

        if test_type == "精度分析":
            half_width_abs = st.number_input(
                "配对差值均值 CI 半宽 w（与均值同单位）",
                value=0.5,
                step=0.05,
                min_value=0.001,
                help="双侧 (1−α) 置信区间对**配对差值均值**的半宽须 ≤ w；按配对差值的 t 半宽法迭代。",
            )
            return {
                "mean_diff": 0.0,
                "sd_diff": sd_diff,
                "margin_percent": 0.0,
                "half_width_abs": half_width_abs,
                "precision_mode": True,
                "paired_precision": True,
            }

        mean_diff = st.number_input("配对差值均值（试验−对照）", value=1.5, step=0.1)
        margin_percent = (
            0.0
            if test_type == "差异性检验"
            else st.number_input("界值 Δ（绝对均值差）", value=1.0, step=0.1)
        )
        paired_method = "t"
        if test_type == "差异性检验":
            _pl = st.selectbox(
                "计算方法（差异性）",
                ["配对 t 检验", "Wilcoxon 符号秩（ARE 校正）"],
                index=0,
                help="优效/非劣/等效固定为配对 t；等效性按 TOST 思路；差异性可在配对 t 与符号秩间选择。",
            )
            paired_method = "wilcoxon" if _pl.startswith("Wilcoxon") else "t"
        return {
            "mean_diff": mean_diff,
            "sd_diff": sd_diff,
            "margin_percent": margin_percent,
            "paired_method": paired_method,
            "precision_mode": False,
            "paired_precision": False,
        }
    if group_type == "三组设计 (Three Arm)":
        if test_type != "差异性检验":
            raise ValueError("三组连续变量当前仅支持差异性检验。")
        st.caption(
            "三组独立样本、假定方差齐性；SD 为共同组内标准差。组1→组3 建议按自然顺序（如安慰剂→低剂量→高剂量）填写。"
        )
        c_a1, c_a2, c_a3 = st.columns(3)
        with c_a1:
            mean1 = st.number_input("组1 预期均值", value=8.0, step=0.1, key="three_m1")
        with c_a2:
            mean2 = st.number_input("组2 预期均值", value=9.5, step=0.1, key="three_m2")
        with c_a3:
            mean3 = st.number_input("组3 预期均值", value=11.0, step=0.1, key="three_m3")
        sd_within = st.number_input(
            "假定共同标准差 SD（组内）",
            value=3.0,
            step=0.1,
            min_value=1e-6,
            key="three_sd",
        )
        _tm = st.selectbox(
            "差异性检验方法",
            [
                "单因素 ANOVA（F 检验）",
                "Kruskal-Wallis 检验（ARE 校正）",
                "有序线性趋势（ANOVA 线性对比）",
                "Dunnett（多试验组 vs 对照）",
            ],
            index=0,
            help="Dunnett：两试验组分别与对照比较；样本量按 Bonferroni 保守估计。",
            key="three_method_sel",
        )
        _tmap = {
            "单因素 ANOVA（F 检验）": "anova_f",
            "Kruskal-Wallis 检验（ARE 校正）": "kw",
            "有序线性趋势（ANOVA 线性对比）": "linear_trend",
            "Dunnett（多试验组 vs 对照）": "dunnett",
        }
        three_arm_method = _tmap[_tm]
        control_index = 0
        if three_arm_method == "dunnett":
            _cl = st.selectbox(
                "对照组（Dunnett）",
                ["组1", "组2", "组3"],
                index=0,
                help="另两组为试验组，分别与对照做比较。",
                key="three_ctrl",
            )
            control_index = ["组1", "组2", "组3"].index(_cl)
        return {
            "mean1": mean1,
            "mean2": mean2,
            "mean3": mean3,
            "sd_within": sd_within,
            "three_arm_method": three_arm_method,
            "control_index": control_index,
            "precision_mode": False,
        }
    raise ValueError(f"不支持的连续变量分组: {group_type}")


def render(col_input, col_result, group_type: str, test_type: str) -> None:
    variable_type = "连续变量 (均值)"

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
    # 左侧：表单（所有输入控件 + 提交按钮）
    # ------------------------------------------------------------------
    with col_input:
        form_key = f"continuous_form_{group_type}_{test_type}"
        with st.form(key=form_key):
            st.subheader("临床预期", divider="gray")
            _c = render_clinical_inputs(group_type, test_type)
            mean_t = _c.get("mean_t")
            mean_c = _c.get("mean_c")
            sd_t = _c.get("sd_t")
            sd_c = _c.get("sd_c")
            mean_ref = _c.get("mean_ref")
            sd_single = _c.get("sd_single")
            mean_diff = _c.get("mean_diff")
            sd_diff = _c.get("sd_diff")
            margin_percent = _c.get("margin_percent", 0.0)
            half_width_abs = _c.get("half_width_abs")
            precision_mode = _c.get("precision_mode", False)
            two_arm_precision = _c.get("two_arm_precision", False)
            paired_precision = _c.get("paired_precision", False)
            paired_method = _c.get("paired_method", "t")
            mean1 = _c.get("mean1")
            mean2 = _c.get("mean2")
            mean3 = _c.get("mean3")
            sd_within = _c.get("sd_within")
            three_arm_method = _c.get("three_arm_method")
            control_index = _c.get("control_index", 0)

            # 统计学设计参数
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
    # 计算
    # ------------------------------------------------------------------
    alloc_ratio_val = eval(alloc_ratio_str.replace(":", "/"))
    continuous_calculator = ContinuousCalculator()
    error_message = None
    result = None
    try:
        if group_type == "独立两组 (Parallel)" and test_type == "精度分析" and two_arm_precision:
            result = continuous_calculator.compute_two_arm_precision(
                sd_t=sd_t,
                sd_c=sd_c,
                half_width=half_width_abs,
                alpha=alpha,
                allocation_ratio=alloc_ratio_val,
                dropout_rate=dropout_percent / 100,
            )
        elif group_type == "独立两组 (Parallel)":
            result = continuous_calculator.compute_two_arm(
                mean_t=mean_t,
                mean_c=mean_c,
                sd_t=sd_t,
                sd_c=sd_c,
                alpha=alpha,
                power=power_percent / 100,
                margin=margin_percent,
                allocation_ratio=alloc_ratio_val,
                dropout_rate=dropout_percent / 100,
                test_type=test_type,
                method=_c.get("two_arm_method", "pooled_t"),
            )
        elif group_type == "单臂设计 (Single Arm)" and test_type == "精度分析" and precision_mode:
            result = continuous_calculator.compute_single_arm_precision(
                sd=sd_single,
                half_width=half_width_abs,
                alpha=alpha,
                dropout_rate=dropout_percent / 100,
            )
        elif group_type == "单臂设计 (Single Arm)":
            result = continuous_calculator.compute_single_arm(
                mean_t=mean_t,
                mean_ref=mean_ref,
                sd=sd_single,
                alpha=alpha,
                power=power_percent / 100,
                margin=margin_percent,
                dropout_rate=dropout_percent / 100,
                test_type=test_type,
                method=_c.get("single_arm_method", "t"),
            )
        elif group_type == "配对两组 (Paired)" and test_type == "精度分析" and paired_precision:
            result = continuous_calculator.compute_paired_precision(
                sd_diff=sd_diff,
                half_width=half_width_abs,
                alpha=alpha,
                dropout_rate=dropout_percent / 100,
            )
        elif group_type == "配对两组 (Paired)":
            result = continuous_calculator.compute_paired(
                mean_diff=mean_diff,
                sd_diff=sd_diff,
                alpha=alpha,
                power=power_percent / 100,
                margin=margin_percent,
                dropout_rate=dropout_percent / 100,
                test_type=test_type,
                method=paired_method,
            )
        elif group_type == "三组设计 (Three Arm)":
            if three_arm_method == "anova_f":
                result = continuous_calculator.compute_three_arm_anova(
                    mean1=mean1,
                    mean2=mean2,
                    mean3=mean3,
                    sd_within=sd_within,
                    alpha=alpha,
                    power=power_percent / 100,
                    dropout_rate=dropout_percent / 100,
                )
            elif three_arm_method == "kw":
                result = continuous_calculator.compute_three_arm_kruskal_wallis(
                    mean1=mean1,
                    mean2=mean2,
                    mean3=mean3,
                    sd_within=sd_within,
                    alpha=alpha,
                    power=power_percent / 100,
                    dropout_rate=dropout_percent / 100,
                )
            elif three_arm_method == "linear_trend":
                result = continuous_calculator.compute_three_arm_linear_trend(
                    mean1=mean1,
                    mean2=mean2,
                    mean3=mean3,
                    sd_within=sd_within,
                    alpha=alpha,
                    power=power_percent / 100,
                    dropout_rate=dropout_percent / 100,
                )
            elif three_arm_method == "dunnett":
                result = continuous_calculator.compute_three_arm_dunnett(
                    mean1=mean1,
                    mean2=mean2,
                    mean3=mean3,
                    control_index=int(control_index),
                    sd_within=sd_within,
                    alpha=alpha,
                    power=power_percent / 100,
                    dropout_rate=dropout_percent / 100,
                )
            else:
                raise ValueError(f"未知的三组连续变量方法: {three_arm_method}")
        else:
            raise ValueError("该连续变量分组设计暂未开放。")
    except Exception as e:
        error_message = str(e)

    # ------------------------------------------------------------------
    # 右侧：结果输出
    # ------------------------------------------------------------------
    with col_result:
        if error_message:
            st.error(f"计算冲突: {error_message}")
            return
        if not result:
            st.warning("无计算结果。")
            return

        if group_type == "三组设计 (Three Arm)":
            total_n = result["total_sample_size_with_dropout"]
            theo_total = result["total_sample_size"]
            theo_pg = int(result["n_per_group"])
            npg_drop = int(result["n_per_group_with_dropout"])

            # 1. Morandi KPI 卡片（三组并列）
            st.markdown(
                f"""
<div class="result-kpi-card">
    <div style="font-size:16px;font-weight:500;text-align:center;opacity:0.95;letter-spacing:0.5px;">
        目标招募样本量
    </div>
    <div class="group-split">
        <div class="group-box">
            <div class="group-num">{npg_drop}</div>
            <div class="group-name">第 1 组</div>
        </div>
        <div class="group-box">
            <div class="group-num">{npg_drop}</div>
            <div class="group-name">第 2 组</div>
        </div>
        <div class="group-box">
            <div class="group-num">{npg_drop}</div>
            <div class="group-name">第 3 组</div>
        </div>
    </div>
    <div class="total-box">
        总计招募：{total_n} 例
        （每组 {npg_drop} 例，已计入 {dropout_percent}% 脱落率）
    </div>
</div>
""",
                unsafe_allow_html=True,
            )

            st.markdown("<br>", unsafe_allow_html=True)
            st.markdown('<div class="section-title">样本量计算方法</div>', unsafe_allow_html=True)

            # 2. 极简版叙述
            ctrl_name = ["组1", "组2", "组3"][int(control_index or 0)]
            if three_arm_method == "dunnett":
                cn_text = (
                    f"本研究为三组连续变量终点，计划采用 Dunnett 型分析（{ctrl_name} 为对照，另两组为试验组分别与对照比较）。"
                    f"假定三组预期均值分别约为 μ₁={mean1:.2f}、μ₂={mean2:.2f}、μ₃={mean3:.2f}，共同组内标准差 SD≈{sd_within:.2f}。"
                    f"采用{alpha_side_label}，显著性水平 α = {alpha}，检验效能 1-β = {power_percent:.0f}%，三组等额分配。"
                    f"在上述参数条件下，依据 {result.get('method', '')} 估算，"
                    f"理论所需总样本量为 {theo_total} 例（每组 {theo_pg} 例）。"
                    f"考虑约 {dropout_percent:.0f}% 脱落率，计划总入组 {total_n} 例。"
                )
            else:
                cn_text = (
                    f"本研究为三组独立样本连续变量终点差异性分析。"
                    f"假定三组预期均值分别约为 μ₁={mean1:.2f}、μ₂={mean2:.2f}、μ₃={mean3:.2f}，"
                    f"共同组内标准差 SD≈{sd_within:.2f}。"
                    f"采用{alpha_side_label}，显著性水平 α = {alpha}，检验效能 1-β = {power_percent:.0f}%，三组等额分配。"
                    f"在上述参数条件下，依据 {result.get('method', '')} 估算，"
                    f"理论所需总样本量为 {theo_total} 例（每组 {theo_pg} 例）。"
                    f"考虑约 {dropout_percent:.0f}% 脱落率，计划总入组 {total_n} 例。"
                )

            st.markdown(
                f"""
<div class="statement-box">
    <p>{cn_text}</p>
</div>
""",
                unsafe_allow_html=True,
            )
            return

        total_n = result["total_sample_size_with_dropout"]
        n_treat = result.get("n_treatment_with_dropout", 0)
        n_ctrl = result.get("n_control_with_dropout", 0)
        theo_total = result["total_sample_size"]

        if group_type in ("单臂设计 (Single Arm)", "配对两组 (Paired)"):
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

            # 2. 极简版叙述
            if group_type == "配对两组 (Paired)" and test_type == "精度分析":
                cn_text = (
                    f"本研究为配对两组连续变量精度分析。"
                    f"假定配对差值标准差 SD_diff≈{sd_diff:.4f}。"
                    f"要求 {(1 - alpha) * 100:.0f}% 置信区间下配对差值均值的半宽不超过 {half_width_abs:.4f}。"
                    f"在上述参数条件下，依据 {result.get('method', 't 分布半宽法')} 计算，"
                    f"理论最小对子数为 {theo_total}。"
                    f"考虑约 {dropout_percent:.0f}% 脱落率，计划入组 {total_n} 例。"
                )
            elif group_type == "配对两组 (Paired)":
                margin_txt = f"、界值 Δ = {margin_percent:.2f}" if test_type != "差异性检验" else ""
                cn_text = (
                    f"本研究为配对两组连续变量{test_type[:2]}设计。"
                    f"假定配对差值均值为 {mean_diff:.2f}，"
                    f"配对差值标准差 SD_diff = {sd_diff:.4f}{margin_txt}。"
                    f"采用{alpha_side_label}，显著性水平 α = {alpha}，检验效能 1-β = {power_percent:.0f}%。"
                    f"在上述参数条件下，依据 {result.get('method', '配对 t 检验')} 计算，"
                    f"理论所需对子数为 {theo_total} 例。"
                    f"考虑约 {dropout_percent:.0f}% 脱落率，实际计划入组 {total_n} 例。"
                )
            elif test_type == "精度分析":
                cn_text = (
                    f"本研究为单臂连续变量精度分析。"
                    f"假定总体标准差 SD≈{sd_single:.2f}，"
                    f"要求 {(1 - alpha) * 100:.0f}% 置信区间下半宽不超过 {half_width_abs:.4f}。"
                    f"在上述参数条件下，依据 {result.get('method', 't 分布半宽法')} 计算，"
                    f"理论最小样本量为 {theo_total} 例。"
                    f"考虑约 {dropout_percent:.0f}% 脱落率，计划入组 {total_n} 例。"
                )
            else:
                margin_txt = f"、界值 Δ = {margin_percent:.2f}" if test_type != "差异性检验" else ""
                cn_text = (
                    f"本研究为单臂连续变量{test_type[:2]}设计。"
                    f"假定试验组均值为 {mean_t:.2f}，参考均值为 {mean_ref:.2f}，"
                    f"标准差 SD 为 {sd_single:.2f}{margin_txt}。"
                    f"采用{alpha_side_label}，显著性水平 α = {alpha}，检验效能 1-β = {power_percent:.0f}%。"
                    f"在上述参数条件下，依据 {result.get('method', '单样本 t 检验')} 计算，"
                    f"理论所需样本量为 {theo_total} 例。"
                    f"考虑约 {dropout_percent:.0f}% 脱落率，实际计划入组 {total_n} 例。"
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

            # 2. 极简版叙述
            if test_type == "精度分析" and two_arm_precision:
                cn_text = (
                    f"本研究为独立两组连续变量均值差精度分析。"
                    f"假定试验组 SD≈{sd_t:.2f}、对照组 SD≈{sd_c:.2f}，"
                    f"要求 {(1 - alpha) * 100:.0f}% 置信区间下均值差的半宽不超过 {half_width_abs:.4f}。"
                    f"在上述参数条件下，依据 {result.get('method', 't 分布半宽法')} 计算，"
                    f"理论最小总样本量为 {theo_total} 例。"
                    f"考虑约 {dropout_percent:.0f}% 脱落率，计划总入组 {total_n} 例。"
                )
            else:
                margin_txt = f"、界值 Δ = {margin_percent:.2f}" if test_type != "差异性检验" else ""
                cn_text = (
                    f"本研究为连续变量终点的独立两组{test_type[:2]}设计。"
                    f"假定试验组均值为 {mean_t:.2f}，对照组均值为 {mean_c:.2f}，"
                    f"试验组标准差 {sd_t:.2f}、对照组标准差 {sd_c:.2f}{margin_txt}。"
                    f"采用{alpha_side_label}，显著性水平 α = {alpha}，检验效能 1-β = {power_percent:.0f}%。"
                    f"在上述参数条件下，依据 {result.get('method', '两独立样本 t')} 计算，"
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
