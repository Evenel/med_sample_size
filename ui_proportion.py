"""分类变量 (率)：临床预期 + 统计学设计 + 调度 + 结果。"""

from __future__ import annotations

import math
import streamlit as st

import ui_common
from calc_proportion import ProportionCalculator, ProportionInput


def render_clinical_inputs(group_type: str, test_type: str, simon_design_type: str | None) -> dict:
    """仅渲染「率」相关临床预期。"""
    if group_type == "Simon 二阶段":
        c_p1, c_p2 = st.columns(2)
        with c_p1:
            p_treatment_percent = st.number_input("目标有效率 P1 (%)", value=40.0, step=0.1)
        with c_p2:
            p_control_percent = st.number_input("无效界值 P0 (%)", value=20.0, step=0.1)
        return {
            "p_treatment_percent": p_treatment_percent,
            "p_control_percent": p_control_percent,
            "margin_percent": 0.0,
        }

    if group_type == "配对两组 (Paired)":
        st.caption(
            "配对二分类：π₁₀=试验有效且对照无效，π₀₁=试验无效且对照有效；边际率差 pT−pC=π₁₀−π₀₁。"
        )
        input_mode = st.radio(
            "参数输入方式",
            ["直接输入 π₁₀、π₀₁（%）", "由边际率与不一致比例 pD 推算"],
            horizontal=True,
        )
        if input_mode.startswith("直接"):
            p10_pct = st.number_input("π₁₀ (%)（试验有效、对照无效）", value=12.0, step=0.1, min_value=0.01)
            p01_pct = st.number_input("π₀₁ (%)（试验无效、对照有效）", value=8.0, step=0.1, min_value=0.01)
        else:
            pt_pct = st.number_input("试验组边际有效率 pT (%)", value=55.0, step=0.1)
            pc_pct = st.number_input("对照组边际有效率 pC (%)", value=50.0, step=0.1)
            pd_pct = st.number_input(
                "不一致对合计比例 pD (%)",
                value=25.0,
                step=0.1,
                min_value=0.1,
                help="须大于 |pT−pC|，否则 π₁₀ 或 π₀₁ 会为非正。",
            )
            psi = (pt_pct - pc_pct) / 100.0
            pdv = pd_pct / 100.0
            p10_pct = 100.0 * (pdv + psi) / 2.0
            p01_pct = 100.0 * (pdv - psi) / 2.0
            st.caption(
                f"推算：π₁₀≈{p10_pct:.2f}%，π₀₁≈{p01_pct:.2f}%（若不合理请调整 pD 或边际率）。"
            )
        direction = "higher"
        if test_type == "差异性检验":
            margin_percent = 0.0
        else:
            margin_label = "非劣界值 Δ" if test_type == "非劣效性检验" else "等效/优效界值 Δ"
            margin_percent = st.number_input(
                f"率差界值 {margin_label} (百分点，ψ=π₁₀−π₀₁)",
                value=5.0,
                step=0.1,
            )
        with st.expander("配对两组率算法说明（严谨版）", expanded=False):
            st.markdown(
                """
| 检验目标 | 场景 | 推荐算法 | 说明 |
|---------|------|----------|------|
| 差异性检验 | 正常样本量 | **McNemar 检验** | 基于 discordant pairs 的正态近似 |
|  | 小样本 / discordant 很少 | **Exact McNemar** | 条件二项精确检验 |
| 优效 / 非劣效 | 常规样本量 | **Nam score** | 配对二分类优/非劣常用标准方法 |
|  | 快速近似 | **paired Wald-type** | 仅粗略估算，不建议作为主方法 |
|  | 极小样本 | **exact conditional paired** | 补充/敏感性分析，非默认主方法 |
| 等效性（TOST） | 常规样本量 | **Nam score TOST** | 双单侧 Nam 型正态近似 |
"""
            )
        return {
            "p10_pct": p10_pct,
            "p01_pct": p01_pct,
            "margin_percent": margin_percent,
            "direction": direction,
        }

    if group_type == "三组设计 (Three Arm)":
        c3a, c3b, c3c = st.columns(3)
        with c3a:
            p1_pct = st.number_input("第 1 组有效率 π₁ (%)", value=50.0, step=0.1)
        with c3b:
            p2_pct = st.number_input("第 2 组有效率 π₂ (%)", value=40.0, step=0.1)
        with c3c:
            p3_pct = st.number_input("第 3 组有效率 π₃ (%)", value=35.0, step=0.1)
        margin_percent = 0.0
        direction = "higher"
        if test_type == "差异性检验（率齐性）":
            st.caption(
                "三组独立二分类：率齐性检验 H₀：π₁=π₂=π₃（任意组间存在差异即可能拒绝）。等额分配；"
                "常规用 Pearson χ² 非中心；若期望频数不足或率极偏，可用蒙特卡洛。"
            )
            with st.expander("三组率齐性检验（算法说明）", expanded=False):
                st.markdown(
                    """
| 场景 | 推荐 | 说明 |
|------|------|------|
| 期望频数较充分 | **Pearson χ² 齐性**（非中心 χ² 样本量） | 与 2×2 表 Pearson χ² 同一渐近族 |
| 小样本 / 极端率 / 需模拟复核 | **蒙特卡洛**（Pearson 统计量 + 渐近 p） | 模拟 3×2 表估计功效，注册试验建议统计师复核 |

* **自动**：三组率均在 [20%,80%] 且按 Pearson 解得的每组最小期望频数 ≥5 → Pearson；否则 → 蒙特卡洛。*
"""
                )
        else:
            st.caption(
                "三组有序（如安慰剂→低剂量→高剂量）：**Cochran–Armitage** 线性趋势检验（双侧）。"
                "分数须严格递增且与临床顺序一致；样本量由蒙特卡洛估计（渐近 Z 检验 p 值）。"
            )
            cs1, cs2, cs3 = st.columns(3)
            with cs1:
                ca_t1 = st.number_input("第 1 组趋势分数", value=0.0, step=0.5)
            with cs2:
                ca_t2 = st.number_input("第 2 组趋势分数", value=1.0, step=0.5)
            with cs3:
                ca_t3 = st.number_input("第 3 组趋势分数", value=2.0, step=0.5)
            with st.expander("Cochran–Armitage 趋势检验说明", expanded=False):
                st.markdown(
                    """
- **适用**：组别有**自然顺序**（剂量、暴露等级等），且主要关心**线性趋势**。
- **不适用**：三组无顺序或效应可能呈非单调（如中剂量最优）——请用侧栏「率齐性」。
- **分数**：默认等距 0,1,2；也可改为与剂量成比例的数值，须**严格递增**。
- **样本量**：蒙特卡洛模拟二项数据 + 双侧渐近 p 值；大样本可与统计师用正态近似交叉验证。
"""
                )
        return {
            "p1_pct": p1_pct,
            "p2_pct": p2_pct,
            "p3_pct": p3_pct,
            "margin_percent": margin_percent,
            "direction": direction,
            "ca_t1": ca_t1 if test_type != "差异性检验（率齐性）" else 0.0,
            "ca_t2": ca_t2 if test_type != "差异性检验（率齐性）" else 1.0,
            "ca_t3": ca_t3 if test_type != "差异性检验（率齐性）" else 2.0,
        }

    if test_type == "精度分析":
        if group_type == "独立两组 (Parallel)":
            c_p1, c_p2 = st.columns(2)
            with c_p1:
                p_treatment_percent = st.number_input("试验组预期有效率 (%)", value=50.0, step=0.1)
            with c_p2:
                p_control_percent = st.number_input("对照组预期有效率 (%)", value=40.0, step=0.1)
            half_width_percent = st.number_input(
                "率差容许误差半宽 w (%)",
                value=10.0,
                step=0.5,
                help="指两组率差 (pT−pC) 的置信区间半宽，例如 10 表示半宽 0.10",
            )
        else:
            p_treatment_percent = st.number_input("预期有效率 P (%)", value=40.0, step=0.1)
            p_control_percent = 0.0
            half_width_percent = st.number_input(
                "容许误差半宽 w (%)",
                value=10.0,
                step=0.5,
                help="例如 10 表示 ±10%，即 95%CI 半宽不超过 0.10",
            )
        margin_percent = 0.0
        direction = "higher"
        return {
            "p_treatment_percent": p_treatment_percent,
            "p_control_percent": p_control_percent,
            "half_width_percent": half_width_percent,
            "margin_percent": margin_percent,
            "direction": direction,
        }

    c_p1, c_p2 = st.columns(2)
    with c_p1:
        p_treatment_percent = st.number_input("试验组预期有效率 (%)", value=60.0, step=0.1)
    label_ctrl = "对照组预期有效率 (%)" if group_type == "独立两组 (Parallel)" else "历史目标有效率 P0 (%)"
    with c_p2:
        p_control_percent = st.number_input(label_ctrl, value=40.0, step=0.1)

    if group_type == "独立两组 (Parallel)" and test_type == "优效性检验":
        direction_label = st.selectbox(
            "指标方向（优效）",
            ["高值更优（有效率更高）", "低值更优（事件率更低）"],
            index=0,
        )
        direction = "lower" if direction_label.startswith("低值更优") else "higher"
    else:
        direction = "higher"

    if test_type == "差异性检验":
        margin_percent = 0.0
    else:
        margin_label = "非劣界值 Δ" if test_type == "非劣效性检验" else "等效/优效界值 Δ"
        margin_percent = st.number_input(f"{margin_label} (百分点)", value=5.0, step=0.1)

    if group_type == "独立两组 (Parallel)":
        st.caption(
            "双臂率：差异性在期望频数≥5时常用合并方差 Z（与 Pearson χ² 等价）；"
            "极端率/小样本倾向 Fisher 或反正弦。优/非劣/等效以 FM Score 为药政常用基准，"
            "等效性对应 FM TOST；可选 Wald 简化或 Chan 类无条件（近似）。精度分析见 Wald / Newcombe。"
        )
        with st.expander("双臂率算法说明（ICH / 临床研究常用）", expanded=False):
            st.markdown(
                """
| 检验目标 | 场景 | 推荐算法 | 备选/补充 |
|---------|------|----------|-----------|
| 差异性（Δ=0） | 期望频数 ≥5（Cochran） | **Pooled Z（合并方差）** | 与 Pearson χ² 等价 |
|  | 极端率 / 小样本 | **Fisher 精确** | 控制 I 类错误 |
| 优效 / 非劣 / 等效 | 期望频数较充分 | **FM Score（边界投影）** | 药政申报常用 |
|  | 极端率 / 极小样本 | **Chan 无条件精确**（本工具为 FM 近似） | StatXact 等复核 |
|  | 大样本简化 | **Wald（非合并边界方差）** | 与 FM 接近 |
| 等效性（TOST） | 期望频数较充分 | **FM TOST** | 两个单侧 FM |
| 精度分析（CI） | 率差置信区间 | **Wald** / **Newcombe-Wilson** | 极端率更稳 |

* **自动（差异性）**：按合并方差估计样本量后，若四格最小期望频数 ≥5（Cochran）→ Pooled Z；否则 → Fisher 精确。罕见病若 n 足够大仍可走 Pooled Z；极端率且期望频数不足时请手动选 Fisher 或反正弦。*
"""
            )

    return {
        "p_treatment_percent": p_treatment_percent,
        "p_control_percent": p_control_percent,
        "margin_percent": margin_percent,
        "direction": direction,
    }


def render(col_input, col_result, group_type: str, test_type: str, simon_design_type: str | None) -> None:
    variable_type = "分类变量 (率)"
    with col_input:
        st.markdown('<div class="section-title">临床预期</div>', unsafe_allow_html=True)
        clinical = render_clinical_inputs(group_type, test_type, simon_design_type)
        p_treatment_percent = clinical.get("p_treatment_percent")
        p_control_percent = clinical.get("p_control_percent")
        p10_pct = clinical.get("p10_pct")
        p01_pct = clinical.get("p01_pct")
        p1_pct = clinical.get("p1_pct")
        p2_pct = clinical.get("p2_pct")
        p3_pct = clinical.get("p3_pct")
        ca_t1 = clinical.get("ca_t1", 0.0)
        ca_t2 = clinical.get("ca_t2", 1.0)
        ca_t3 = clinical.get("ca_t3", 2.0)
        margin_percent = clinical.get("margin_percent", 0.0)
        direction = clinical.get("direction", "higher")
        half_width_percent = clinical.get("half_width_percent")

        stats = ui_common.render_statistics_design(variable_type, group_type, test_type)
        alpha = stats["alpha"]
        power_percent = stats["power_percent"]
        dropout_percent = stats["dropout_percent"]
        alloc_ratio_str = stats["alloc_ratio_str"]
        precision_method = stats["precision_method"]
        single_arm_method = stats["single_arm_method"]
        two_arm_diff_method = stats["two_arm_diff_method"]
        two_arm_margin_method = stats["two_arm_margin_method"]
        paired_rate_method = stats["paired_rate_method"]
        three_arm_diff_method = stats["three_arm_diff_method"]
        alpha_side_label = stats["alpha_side_label"]
        if half_width_percent is None:
            half_width_percent = stats.get("half_width_percent", 10.0)

    alloc_ratio_val = eval(alloc_ratio_str.replace(":", "/"))
    proportion_calculator = ProportionCalculator()
    error_message = None
    result = None
    simon_res = None
    try:
        if group_type == "Simon 二阶段":
            simon_res = proportion_calculator.compute_simon_two_stage(
                p0=p_control_percent / 100,
                p1=p_treatment_percent / 100,
                alpha=alpha,
                power=power_percent / 100,
                dropout_rate=dropout_percent / 100,
                design_type=simon_design_type or "optimal",
            )
        elif (
            group_type == "独立两组 (Parallel)"
            and test_type == "精度分析"
        ):
            pm = "newcombe_wilson" if precision_method.startswith("Newcombe") else "wald"
            prec = proportion_calculator.compute_two_arm_precision(
                p_t=p_treatment_percent / 100,
                p_c=p_control_percent / 100,
                r=alloc_ratio_val,
                ci_level=1 - alpha,
                half_width=half_width_percent / 100,
                method=pm,
            )
            nt = int(prec["n_treatment"])
            nc = int(prec["n_control"])
            ntd = int(math.ceil(nt / (1 - dropout_percent / 100)))
            ncd = int(math.ceil(nc / (1 - dropout_percent / 100)))
            result = {
                "n_treatment": nt,
                "n_control": nc,
                "total_sample_size": nt + nc,
                "n_treatment_with_dropout": ntd,
                "n_control_with_dropout": ncd,
                "total_sample_size_with_dropout": ntd + ncd,
                "method": prec["method"],
                "precision": prec,
            }
        elif group_type == "单臂设计 (Single Arm)" and test_type == "精度分析":
            prec_method = "clopper_pearson" if precision_method.startswith("确切") else "wald"
            prec = proportion_calculator.compute_single_arm_precision(
                p=p_treatment_percent / 100,
                ci_level=1 - alpha,
                half_width=half_width_percent / 100,
                method=prec_method,
            )
            n = int(prec["n"])
            n_drop = int(math.ceil(n / (1 - dropout_percent / 100)))
            result = {
                "n_treatment": n,
                "n_control": 0,
                "total_sample_size": n,
                "n_treatment_with_dropout": n_drop,
                "n_control_with_dropout": 0,
                "total_sample_size_with_dropout": n_drop,
                "method": prec["method"],
                "precision": prec,
            }
        elif group_type == "配对两组 (Paired)":
            result = proportion_calculator.compute_paired_proportion(
                p_discordant_treat=p10_pct / 100.0,
                p_discordant_ctrl=p01_pct / 100.0,
                alpha=alpha,
                power=power_percent / 100.0,
                dropout_rate=dropout_percent / 100.0,
                test_type=test_type,
                margin=margin_percent / 100.0 if test_type != "差异性检验" else 0.0,
                method=paired_rate_method,
            )
        elif group_type == "三组设计 (Three Arm)":
            if test_type == "差异性检验（率齐性）":
                result = proportion_calculator.compute_three_arm_proportion(
                    p1_pct / 100.0,
                    p2_pct / 100.0,
                    p3_pct / 100.0,
                    alpha=alpha,
                    power=power_percent / 100.0,
                    dropout_rate=dropout_percent / 100.0,
                    method=three_arm_diff_method,
                )
            else:
                result = proportion_calculator.compute_three_arm_cochran_armitage(
                    p1_pct / 100.0,
                    p2_pct / 100.0,
                    p3_pct / 100.0,
                    alpha=alpha,
                    power=power_percent / 100.0,
                    dropout_rate=dropout_percent / 100.0,
                    scores=(ca_t1, ca_t2, ca_t3),
                )
        else:
            params = ProportionInput(
                p_treatment=p_treatment_percent / 100,
                p_control=p_control_percent / 100,
                alpha=alpha,
                power=power_percent / 100,
                test_type=test_type,
                allocation_ratio=alloc_ratio_val,
                dropout_rate=dropout_percent / 100,
                margin=margin_percent / 100,
                direction=direction if group_type == "独立两组 (Parallel)" else "higher",
                two_arm_diff_method=two_arm_diff_method,
                two_arm_margin_method=two_arm_margin_method,
            )
            if group_type == "独立两组 (Parallel)":
                result = proportion_calculator.compute_two_arm(params)
            else:
                result = proportion_calculator.compute_single_arm(params, method=single_arm_method)
    except Exception as e:
        error_message = str(e)

    with col_result:
        if error_message:
            st.error(f"计算冲突: {error_message}")
            return

        if group_type == "Simon 二阶段" and simon_res:
            st.markdown(
                f"""
                <div class="result-kpi-card">
                    <div style="font-size: 16px; font-weight: 500; text-align: center; opacity: 0.95;">Simon {simon_res['design_type'].capitalize()} 目标招募样本量</div>
                    <div class="group-split">
                        <div class="group-box">
                            <div class="group-num">{simon_res['n1']}</div>
                            <div class="group-name">第一阶段</div>
                        </div>
                        <div class="group-box">
                            <div class="group-num">{simon_res['n2']}</div>
                            <div class="group-name">第二阶段</div>
                        </div>
                    </div>
                    <div class="total-box">
                        理论总样本量：{simon_res['n']} 例；实际计划入组：{simon_res['total_sample_size_with_dropout']} 例
                        <span style="font-weight: 400; font-size: 15px; opacity: 0.9;">(已计入 {dropout_percent}% 脱落率)</span>
                    </div>
                </div>
                """,
                unsafe_allow_html=True,
            )
            st.caption(
                f"设计参数：r1={simon_res['r1']}，r={simon_res['r']}，"
                f"PET(p0)={simon_res['pet0']*100:.1f}%，EN(p0)={simon_res['en0']:.1f}，"
                f"实际α={simon_res['actual_alpha']:.4f}，实际Power={simon_res['actual_power']:.4f}"
            )
            st.caption(
                f"对照设计：Optimal(n={simon_res['optimal_design']['n']}, EN0={simon_res['optimal_design']['en0']:.1f})；"
                f"Minimax(n={simon_res['minimax_design']['n']}, EN0={simon_res['minimax_design']['en0']:.1f})"
            )
            st.markdown('<div class="section-title" style="margin-top:0;">样本量计算方法</div>', unsafe_allow_html=True)
            cn_text = (
                f"本研究采用 Simon 二阶段 {simon_res['design_type']} 设计。设定药物无效界值 P0 为 {p_control_percent}%，"
                f"预期目标有效率 P1 为 {p_treatment_percent}%。在单侧显著性水平 α = {alpha} 且检验效能 1-β = {power_percent}% 的条件下，"
                f"第一阶段计划入组 {simon_res['n1']} 例受试者。若在第一阶段观察到 ≤ {simon_res['r1']} 例有效，则因无效提前终止试验"
                f"（早期停止概率 PET = {simon_res['pet0']*100:.1f}%）；若观察到 > {simon_res['r1']} 例有效，则进入第二阶段，"
                f"继续入组 {simon_res['n2']} 例受试者。<br><br>"
                f"最终，若在总计理论 {simon_res['n']} 例受试者中观察到 > {simon_res['r']} 例有效，则拒绝原假设，"
                f"认为该药物具有临床疗效。考虑 {dropout_percent}% 的脱落率，实际计划入组最大样本量为 {simon_res['total_sample_size_with_dropout']} 例。"
            )
            st.markdown(f"<div class='statement-box'><p>{cn_text}</p></div>", unsafe_allow_html=True)
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

        if group_type == "三组设计 (Three Arm)":
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
            if result.get("method"):
                st.caption(result["method"])
            if result.get("rate_note"):
                st.caption(result["rate_note"])
            if result.get("min_expected_cell_approx") is not None:
                st.caption(
                    f"按假设率估计的齐性检验最小期望频数（每组阳性/阴性中较小者）≈ {result['min_expected_cell_approx']:.2f}"
                )
            if result.get("actual_power_mc") is not None:
                st.caption(
                    f"蒙特卡洛估计功效 ≈ {result['actual_power_mc'] * 100:.2f}%（{result.get('mc_replications', '')} 次重复）"
                )
            st.markdown('<div class="section-title" style="margin-top:0;">样本量计算方法</div>', unsafe_allow_html=True)
            if test_type == "差异性检验（率齐性）":
                cn_three_1 = (
                    f"本研究为三组独立样本二分类终点**率齐性**分析。假定三组有效率分别约为 π₁={p1_pct:.1f}%、"
                    f"π₂={p2_pct:.1f}%、π₃={p3_pct:.1f}%。"
                    f"采用双侧检验，显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%，三组等额分配。"
                )
            else:
                sc = result.get("ca_scores") or (0.0, 1.0, 2.0)
                cn_three_1 = (
                    f"本研究为三组有序二分类终点的**线性趋势**分析（Cochran–Armitage，双侧）。"
                    f"假定三组有效率分别约为 π₁={p1_pct:.1f}%、π₂={p2_pct:.1f}%、π₃={p3_pct:.1f}%，"
                    f"对应趋势分数为 ({sc[0]:g}, {sc[1]:g}, {sc[2]:g})。"
                    f"采用双侧检验，显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%，三组等额分配。"
                )
            cn_three_2 = (
                f"样本量采用【{result.get('method', '')}】估算。"
                f"理论所需总样本量为 {theo_total} 例（每组 {theo_pg} 例）；考虑约 {dropout_percent:.1f}% 脱落率，计划总入组 {total_n} 例。"
            )
            st.markdown(
                f"<div class='statement-box'><p>{cn_three_1}</p><p>{cn_three_2}</p></div>",
                unsafe_allow_html=True,
            )
            return

        if group_type in ("单臂设计 (Single Arm)", "配对两组 (Paired)"):
            _kpi_title = "目标配对受试者数" if group_type == "配对两组 (Paired)" else "目标招募样本量"
            _kpi_sub = "每人提供配对二分类（两次测量）" if group_type == "配对两组 (Paired)" else ""
            st.markdown(
                f"""
                <div class="result-kpi-card">
                    <div style="font-size: 16px; font-weight: 500; text-align: center; opacity: 0.95;">{_kpi_title}</div>
                    <div style="text-align: center; margin: 30px 0;">
                        <span style="font-size: 64px; font-weight: 700;">{total_n}</span>
                    </div>
                    <div class="total-box">理论所需：{theo_total} 例 <span style="font-weight: 400; font-size: 15px; opacity: 0.9;">(已计入 {dropout_percent}% 脱落率)</span></div>
                </div>
                """,
                unsafe_allow_html=True,
            )
            if _kpi_sub:
                st.caption(_kpi_sub)
            st.markdown('<div class="section-title" style="margin-top:0;">样本量计算方法</div>', unsafe_allow_html=True)
            if group_type == "配对两组 (Paired)" and test_type != "精度分析":
                margin_txt = (
                    f"，设定率差界值 Δ = {margin_percent:.1f}%"
                    if test_type != "差异性检验"
                    else ""
                )
                cn_text_p1 = (
                    f"本研究为配对两组二分类终点{test_type[:2]}设计。假定不一致配对比例 π₁₀≈{p10_pct:.2f}%、"
                    f"π₀₁≈{p01_pct:.2f}%（ψ=π₁₀−π₀₁）{margin_txt}。"
                    f"采用{alpha_side_label}检验，显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                )
                cn_text_p2 = (
                    f"样本量采用【{result.get('method', '')}】计算。"
                    f"理论所需配对受试者数为 {theo_total} 例；考虑约 {dropout_percent:.1f}% 脱落率，计划入组 {total_n} 例。"
                )
                if isinstance(result, dict) and result.get("rate_note"):
                    st.caption(result["rate_note"])
            elif test_type == "精度分析":
                cn_text_p1 = (
                    f"本研究为单臂率精度分析。假定目标率约为 {p_treatment_percent:.1f}%，"
                    f"设定 {(1 - alpha) * 100:.0f}% 置信区间半宽不超过 {half_width_percent:.1f}%。"
                )
                cn_text_p2 = (
                    f"样本量采用【{result.get('method', '')}】计算。"
                    f"理论最小样本量为 {theo_total} 例；考虑约 {dropout_percent:.1f}% 脱落率，计划入组 {total_n} 例。"
                )
            else:
                margin_txt = f"，设定{test_type[:2]}界值 Δ = {margin_percent:.1f}%" if test_type != "差异性检验" else ""
                method_line = result.get("method", "单臂率：Score 正态近似") if isinstance(result, dict) else "单臂率：Score 正态近似"
                cn_text_p1 = (
                    f"本研究为单臂{test_type[:2]}设计。根据临床预期，假定试验组有效率为 {p_treatment_percent:.1f}%，"
                    f"历史对照目标值为 {p_control_percent:.1f}%{margin_txt}。采用{alpha_side_label}检验，"
                    f"显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                )
                cn_text_p2 = (
                    f"样本量采用【{method_line}】进行估计。在上述参数条件下，理论所需总样本量为 {theo_total} 例。"
                    f"考虑约 {dropout_percent:.1f}% 的脱落率，为保证最终分析所需样本量，实际计划总入组 {total_n} 例。"
                )
            if isinstance(result, dict) and "actual_alpha_at_p0" in result:
                st.caption(
                    f"Exact 校核：实际α(p0)≈{result['actual_alpha_at_p0']:.4f}，"
                    f"实际Power(p1)≈{result.get('actual_power_at_p1', float('nan')):.4f}；"
                    f"搜索起点 n0={result.get('n0_start','-')}，n_start={result.get('n_start','-')}。"
                )
            st.markdown(
                f"<div class='statement-box'><p>{cn_text_p1}</p><p>{cn_text_p2}</p></div>",
                unsafe_allow_html=True,
            )
            return

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
        if group_type == "独立两组 (Parallel)" and result:
            if result.get("method"):
                st.caption(result["method"])
            if result.get("expected_min_cell") is not None:
                st.caption(
                    f"基于假设率估计的四格最小期望频数 ≈ {result['expected_min_cell']:.2f}"
                )
            if result.get("rate_note"):
                st.caption(result["rate_note"])
        st.markdown('<div class="section-title" style="margin-top:0;">样本量计算方法</div>', unsafe_allow_html=True)
        alloc_txt = f"即每组 {theo_treat} 例" if alloc_ratio_val == 1.0 else f"即试验组 {theo_treat} 例，对照组 {theo_ctrl} 例"
        if test_type == "精度分析":
            cn_text_p1 = (
                f"本研究为两臂率差精度分析。假定试验组与对照组有效率分别约为 {p_treatment_percent:.1f}%、{p_control_percent:.1f}%，"
                f"设定 {(1 - alpha) * 100:.0f}% 置信水平下率差置信区间半宽不超过 {half_width_percent:.1f}%。"
            )
            cn_text_p2 = (
                f"样本量采用【{result.get('method', '')}】计算。"
                f"理论最小总样本量为 {theo_total} 例（{alloc_txt}）；考虑约 {dropout_percent:.1f}% 脱落率，计划入组 {total_n} 例。"
            )
        else:
            margin_txt = f"，设定{test_type[:2]}界值 Δ = {margin_percent:.1f}%" if test_type != "差异性检验" else ""
            if result and result.get("rate_zone") == "fisher":
                meth_line = "Fisher 精确检验（两独立样本率，双侧）"
            elif result and result.get("rate_zone") == "extreme" and test_type == "差异性检验":
                meth_line = "反正弦变换（方差稳定化）法（两独立样本率差）"
            elif (
                result
                and result.get("rate_zone") == "extreme_pooled"
                and test_type == "差异性检验"
            ):
                meth_line = "合并方差 Z（与 Pearson χ² 等价；率偏斜但期望频数仍 ≥5）"
            elif (
                result
                and result.get("rate_zone") == "extreme_pooled_sparse"
                and test_type == "差异性检验"
            ):
                meth_line = "合并方差 Z；期望频数不足 5，建议分析阶段 Fisher 或改用自动/反正弦"
            elif result and result.get("rate_zone") == "extreme" and test_type in (
                "优效性检验",
                "非劣效性检验",
                "等效性检验",
            ):
                meth_line = "Farrington–Manning score / FM Score（两独立样本率差，H₀ 边界投影）"
            elif result and result.get("rate_zone") == "wald":
                meth_line = "Wald 型非合并方差（H₀ 边界直接代入）"
            elif test_type == "等效性检验":
                meth_line = "FM TOST（两个单侧 FM Score，正态近似）"
            else:
                meth_line = "正态近似法（两独立样本率差）"
            cn_text_p1 = (
                f"主要终点为二分类结局（有效率）。本研究采用独立两组{test_type[:2]}设计。"
                f"根据既往研究及临床预期，试验组有效率预计为 {p_treatment_percent:.1f}%，对照组为 {p_control_percent:.1f}%{margin_txt}。"
                f"采用{alpha_side_label}检验，显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%，按 {alloc_ratio_str} 比例分配。"
            )
            cn_text_p2 = (
                f"样本量采用【{meth_line}】计算。按上述假设与参数，理论所需总样本量为 {theo_total} 例（{alloc_txt}）。"
                f"考虑约 {dropout_percent:.1f}% 脱落率，为保证达到理论分析样本量，计划最多入组 {total_n} 例，"
                f"其中试验组 {n_treat} 例、对照组 {n_ctrl} 例。"
            )
        st.markdown(
            f"<div class='statement-box'><p>{cn_text_p1}</p><p>{cn_text_p2}</p></div>",
            unsafe_allow_html=True,
        )
