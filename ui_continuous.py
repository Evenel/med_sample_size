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
        st.caption(
            "配对设计以**配对差值**（如试验−对照）为分析单位；SD_diff 为差值的标准差。"
            " 亦可由 SD_T、SD_C 与同一受试者两时点/两臂的**组内相关系数 ρ** 推算："
            " SD_diff = √(SD_T² + SD_C² − 2ρ·SD_T·SD_C)。ρ 常取 0.3–0.7，默认 0.5；若不确定建议做敏感性分析。"
        )
        paired_sd_mode = st.radio(
            "配对差值离散程度",
            ["直接输入 SD_diff", "由 SD_T、SD_C、ρ 推算"],
            horizontal=True,
            index=0,
        )
        sd_t = None
        sd_c = None
        rho = None
        if paired_sd_mode == "直接输入 SD_diff":
            sd_diff = st.number_input(
                "配对差值标准差 SD_diff",
                value=3.0,
                step=0.1,
                min_value=1e-6,
            )
        else:
            c_sd1, c_sd2 = st.columns(2)
            with c_sd1:
                sd_t = st.number_input("试验组标准差 SD_T", value=3.0, step=0.1, min_value=1e-6)
            with c_sd2:
                sd_c = st.number_input("对照组标准差 SD_C", value=3.0, step=0.1, min_value=1e-6)
            rho = st.number_input(
                "组内相关系数 ρ",
                value=0.5,
                step=0.05,
                min_value=-0.999,
                max_value=0.999,
                help="同一受试者上两观测（或配对两臂）的相关系数；ρ 越大则 SD_diff 越小。",
            )
            var_d = sd_t**2 + sd_c**2 - 2.0 * rho * sd_t * sd_c
            if var_d <= 1e-12:
                raise ValueError(
                    "配对差值方差 SD_T² + SD_C² − 2ρ·SD_T·SD_C 须为正，请调整 ρ 或两组 SD。"
                )
            sd_diff = math.sqrt(var_d)
            st.metric("推算得到的配对差值标准差 SD_diff", f"{sd_diff:.4f}")

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
                "sd_t": sd_t,
                "sd_c": sd_c,
                "rho": rho,
                "paired_sd_mode": paired_sd_mode,
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
            "sd_t": sd_t,
            "sd_c": sd_c,
            "rho": rho,
            "paired_sd_mode": paired_sd_mode,
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
    with col_input:
        st.markdown('<div class="section-title">临床预期</div>', unsafe_allow_html=True)
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
        paired_sd_mode = _c.get("paired_sd_mode")
        rho = _c.get("rho")
        sd_t_paired = _c.get("sd_t")
        sd_c_paired = _c.get("sd_c")
        mean1 = _c.get("mean1")
        mean2 = _c.get("mean2")
        mean3 = _c.get("mean3")
        sd_within = _c.get("sd_within")
        three_arm_method = _c.get("three_arm_method")
        control_index = _c.get("control_index", 0)

        stats = ui_common.render_statistics_design(variable_type, group_type, test_type)
        alpha = stats["alpha"]
        power_percent = stats["power_percent"]
        dropout_percent = stats["dropout_percent"]
        alloc_ratio_str = stats["alloc_ratio_str"]
        alpha_side_label = stats["alpha_side_label"]

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
            st.markdown(
                f"""
                <div class="result-kpi-card">
                    <div style="font-size: 16px; font-weight: 500; text-align: center; opacity: 0.95;">目标招募总样本量</div>
                    <div style="text-align: center; margin: 30px 0;">
                        <span style="font-size: 64px; font-weight: 700;">{total_n}</span>
                    </div>
                    <div class="total-box">理论所需：{theo_total} 例（三组等额，每组 {theo_pg} 例）
                        <span style="font-weight: 400; font-size: 15px; opacity: 0.9;">(已计入 {dropout_percent}% 脱落率)</span>
                    </div>
                </div>
                """,
                unsafe_allow_html=True,
            )
            st.caption(f"计划入组：每组约 {npg_drop} 例（含脱落）。")
            if result.get("anova_n_per_group_reference") is not None:
                st.caption(
                    f"参考：同参数下单因素 ANOVA 每组约 {result['anova_n_per_group_reference']} 例；"
                    f"Kruskal-Wallis 按 ARE≈0.955 放大。"
                )
            if result.get("dunnett_min_pairwise_delta") is not None:
                st.caption(
                    f"Dunnett 保守估计所依据的两臂与对照最小绝对均值差 ≈ {result['dunnett_min_pairwise_delta']:.4f}（与 SD 同单位）。"
                )
            st.markdown('<div class="section-title" style="margin-top:0;">样本量计算方法</div>', unsafe_allow_html=True)
            _meth = result.get("method", "") if isinstance(result, dict) else ""
            ctrl_name = ["组1", "组2", "组3"][int(control_index or 0)]
            cn_three_1 = (
                f"本研究为三组独立样本连续变量终点**差异性**分析。假定三组预期均值分别约为 μ₁={mean1:.2f}、μ₂={mean2:.2f}、μ₃={mean3:.2f}，"
                f"共同组内标准差 SD≈{sd_within:.2f}。"
                f"采用{alpha_side_label}检验，显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%，三组等额分配。"
            )
            if three_arm_method == "dunnett":
                cn_three_1 = (
                    f"本研究为三组连续变量终点，计划采用 **Dunnett** 型分析（**{ctrl_name}** 为对照，另两组为试验组分别与对照比较）。"
                    f"假定三组预期均值分别约为 μ₁={mean1:.2f}、μ₂={mean2:.2f}、μ₃={mean3:.2f}，共同组内标准差 SD≈{sd_within:.2f}。"
                    f"采用{alpha_side_label}检验，显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%，三组等额分配。"
                )
            cn_three_2 = (
                f"样本量采用【{_meth}】估算。"
                f"理论所需总样本量为 {theo_total} 例（每组 {theo_pg} 例）；考虑约 {dropout_percent:.1f}% 脱落率，计划总入组 {total_n} 例。"
            )
            st.markdown(
                f"<div class='statement-box'><p>{cn_three_1}</p><p>{cn_three_2}</p></div>",
                unsafe_allow_html=True,
            )
            return

        total_n = result["total_sample_size_with_dropout"]
        n_treat = result.get("n_treatment_with_dropout", 0)
        n_ctrl = result.get("n_control_with_dropout", 0)
        theo_total = result["total_sample_size"]
        theo_treat = int(result.get("n_treatment", theo_total / (1 + 1 / alloc_ratio_val)))
        theo_ctrl = int(result.get("n_control", theo_total - theo_treat))

        if group_type in ("单臂设计 (Single Arm)", "配对两组 (Paired)"):
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
            if (
                group_type == "配对两组 (Paired)"
                and test_type == "精度分析"
                and precision_mode
                and paired_precision
            ):
                cn_text_p1 = (
                    f"本研究为配对两组连续变量**精度分析**。假定配对差值标准差 SD_diff≈{sd_diff:.4f}"
                    + (
                        f"（由 SD_T={sd_t_paired:.2f}、SD_C={sd_c_paired:.2f}、ρ={rho:.2f} 推算）。"
                        if paired_sd_mode == "由 SD_T、SD_C、ρ 推算"
                        and sd_t_paired is not None
                        and sd_c_paired is not None
                        and rho is not None
                        else "。"
                    )
                    + f" 要求 {(1 - alpha) * 100:.0f}% 置信区间下配对差值均值的半宽不超过 {half_width_abs:.4f}（与均值同单位）。"
                )
                _meth = result.get("method", "t 分布半宽法") if isinstance(result, dict) else "t 分布半宽法"
                cn_text_p2 = (
                    f"样本量采用【{_meth}】。理论最小对子数为 {theo_total}；"
                    f"考虑约 {dropout_percent:.1f}% 脱落率，计划入组 {total_n} 例。"
                )
            elif group_type == "配对两组 (Paired)":
                margin_txt = f"，界值 Δ = {margin_percent:.2f}" if test_type != "差异性检验" else ""
                if paired_sd_mode == "由 SD_T、SD_C、ρ 推算" and sd_t_paired is not None and rho is not None:
                    sd_txt = (
                        f"由 SD_T = {sd_t_paired:.2f}、SD_C = {sd_c_paired:.2f}、ρ = {rho:.2f} 推算得"
                        f" SD_diff = {sd_diff:.4f}"
                    )
                else:
                    sd_txt = f"配对差值标准差 SD_diff = {sd_diff:.4f}"
                cn_text_p1 = (
                    f"本研究为配对两组连续变量{test_type[:2]}设计。假定配对差值均值为 {mean_diff:.2f}，"
                    f"{sd_txt}{margin_txt}。采用{alpha_side_label}检验，"
                    f"显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                )
                _meth = (
                    result.get("method", "配对 t 检验")
                    if isinstance(result, dict)
                    else "配对 t 检验"
                )
                cn_text_p2 = (
                    f"样本量采用【{_meth}】。在上述参数条件下，理论所需对子数为 {theo_total} 例。"
                    f"考虑约 {dropout_percent:.1f}% 的脱落率，为保证最终分析所需样本量，实际计划入组 {total_n} 例。"
                )
            elif test_type == "精度分析" and precision_mode:
                cn_text_p1 = (
                    f"本研究为单臂连续变量**精度分析**。假定总体标准差 SD≈{sd_single:.2f}，"
                    f"要求 {(1 - alpha) * 100:.0f}% 置信区间下半宽不超过 {half_width_abs:.4f}（与均值同单位）。"
                )
                _meth = result.get("method", "t 分布半宽法") if isinstance(result, dict) else "t 分布半宽法"
                cn_text_p2 = (
                    f"样本量采用【{_meth}】。理论最小样本量为 {theo_total} 例；"
                    f"考虑约 {dropout_percent:.1f}% 脱落率，计划入组 {total_n} 例。"
                )
            else:
                margin_txt = f"，界值 Δ = {margin_percent:.2f}" if test_type != "差异性检验" else ""
                cn_text_p1 = (
                    f"本研究为单臂连续变量{test_type[:2]}设计。假定试验组均值为 {mean_t:.2f}，"
                    f"参考均值为 {mean_ref:.2f}，标准差 SD 为 {sd_single:.2f}{margin_txt}。"
                    f"采用{alpha_side_label}检验，显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                )
                _meth = (
                    result.get("method", "单样本 t 检验")
                    if isinstance(result, dict)
                    else "单样本 t 检验"
                )
                cn_text_p2 = (
                    f"样本量采用【{_meth}】。在上述参数条件下，理论所需样本量为 {theo_total} 例。"
                    f"考虑约 {dropout_percent:.1f}% 的脱落率，为保证最终分析所需样本量，实际计划入组 {total_n} 例。"
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
            if test_type == "精度分析" and two_arm_precision:
                cn_text_p1 = (
                    f"本研究为独立两组连续变量**均值差精度分析**。假定试验组 SD≈{sd_t:.2f}、对照组 SD≈{sd_c:.2f}，"
                    f"要求 {(1 - alpha) * 100:.0f}% 置信区间下均值差 (μT−μC) 的半宽不超过 {half_width_abs:.4f}（与均值同单位）。"
                )
                _meth = result.get("method", "t 分布半宽法（Welch）") if isinstance(result, dict) else "t 分布半宽法（Welch）"
                cn_text_p2 = (
                    f"样本量采用【{_meth}】。理论最小总样本量为 {theo_total} 例（{alloc_txt}）；"
                    f"考虑约 {dropout_percent:.1f}% 脱落率，计划总入组 {total_n} 例。"
                )
            else:
                margin_txt = f"，界值 Δ = {margin_percent:.2f}" if test_type != "差异性检验" else ""
                cn_text_p1 = (
                    f"本研究为连续变量终点的独立两组{test_type[:2]}设计。假定试验组均值为 {mean_t:.2f}，"
                    f"对照组均值为 {mean_c:.2f}，试验组标准差 {sd_t:.2f}、对照组标准差 {sd_c:.2f}{margin_txt}。"
                    f"采用{alpha_side_label}检验，显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                )
                _meth = (
                    result.get("method", "两独立样本 t（合并方差）")
                    if isinstance(result, dict)
                    else "两独立样本 t（合并方差）"
                )
                cn_text_p2 = (
                    f"样本量采用【{_meth}】。在上述参数条件下，理论总样本量为 {theo_total} 例，"
                    f"{alloc_txt}。考虑约 {dropout_percent:.1f}% 的脱落率，实际计划总入组 {total_n} 例，"
                    f"其中试验组 {n_treat} 例，对照组 {n_ctrl} 例。"
                )
            st.markdown(
                f"<div class='statement-box'><p>{cn_text_p1}</p><p>{cn_text_p2}</p></div>",
                unsafe_allow_html=True,
            )
