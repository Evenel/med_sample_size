"""分类变量 (率)：临床预期 + 统计学设计 + 调度 + 结果输出。

输出文字格式符合临床试验方案（CSP）样本量估算章节的写作规范，
可直接用于方案撰写及伦理申请。
"""

from __future__ import annotations

import math

import streamlit as st

import ui_common
from calc_proportion import ProportionCalculator, ProportionInput


def render_clinical_inputs(group_type: str, test_type: str,
                           simon_design_type: str | None) -> dict:
    """渲染「率」相关临床预期输入区。"""
    # ---- Simon 二阶段 ----
    if group_type == "Simon 二阶段":
        c_p1, c_p2 = st.columns(2)
        with c_p1:
            p_treatment_percent = st.number_input(
                "目标有效率 P1 (%)", value=40.0, step=0.1,
            )
        with c_p2:
            p_control_percent = st.number_input(
                "无效界值 P0 (%)", value=20.0, step=0.1,
            )
        return {
            "p_treatment_percent": p_treatment_percent,
            "p_control_percent": p_control_percent,
            "margin_percent": 0.0,
        }

    # ---- 配对两组 ----
    if group_type == "配对两组 (Paired)":
        st.caption(
            "配对二分类：π₁₀ = 试验有效且对照无效，π₀₁ = 试验无效且对照有效。"
            "边际率差 pT − pC = π₁₀ − π₀₁。"
        )
        p10_pct = st.number_input(
            "π₁₀ (%)（试验有效、对照无效）",
            value=12.0, step=0.1, min_value=0.01,
        )
        p01_pct = st.number_input(
            "π₀₁ (%)（试验无效、对照有效）",
            value=8.0, step=0.1, min_value=0.01,
        )
        direction = "higher"
        margin_percent = 0.0
        if test_type != "差异性检验":
            margin_label = (
                "非劣效界值 Δ" if test_type == "非劣效性检验"
                else "等效/优效界值 Δ"
            )
            margin_percent = st.number_input(
                f"率差界值 {margin_label} (百分点，ψ = π₁₀ − π₀₁)",
                value=5.0, step=0.1,
            )

        return {
            "p10_pct": p10_pct,
            "p01_pct": p01_pct,
            "margin_percent": margin_percent,
            "direction": direction,
        }

    # ---- 三组设计 ----
    if group_type == "三组设计 (Three Arm)":
        c3a, c3b, c3c = st.columns(3)
        with c3a:
            p1_pct = st.number_input(
                "第 1 组有效率 π₁ (%)", value=50.0, step=0.1,
            )
        with c3b:
            p2_pct = st.number_input(
                "第 2 组有效率 π₂ (%)", value=40.0, step=0.1,
            )
        with c3c:
            p3_pct = st.number_input(
                "第 3 组有效率 π₃ (%)", value=35.0, step=0.1,
            )
        margin_percent = 0.0
        direction = "higher"
        ca_t1, ca_t2, ca_t3 = 0.0, 1.0, 2.0

        if test_type == "差异性检验（率齐性）":
            st.caption(
                "三组独立二分类：率齐性检验 H₀: π₁ = π₂ = π₃。"
                "等额分配，采用 Pearson χ² 非中心法。"
            )
        else:
            st.caption(
                "三组有序（如安慰剂 → 低剂量 → 高剂量）："
                "**Cochran–Armitage** 线性趋势检验（双侧）。"
                "分数须严格递增且与临床顺序一致；采用非中心 χ² 法估计样本量。"
            )
            cs1, cs2, cs3 = st.columns(3)
            with cs1:
                ca_t1 = st.number_input(
                    "第 1 组趋势分数", value=0.0, step=0.5,
                )
            with cs2:
                ca_t2 = st.number_input(
                    "第 2 组趋势分数", value=1.0, step=0.5,
                )
            with cs3:
                ca_t3 = st.number_input(
                    "第 3 组趋势分数", value=2.0, step=0.5,
                )
        return {
            "p1_pct": p1_pct, "p2_pct": p2_pct, "p3_pct": p3_pct,
            "margin_percent": margin_percent, "direction": direction,
            "ca_t1": ca_t1, "ca_t2": ca_t2, "ca_t3": ca_t3,
        }

    # ---- 精度分析 ----
    if test_type == "精度分析":
        if group_type == "独立两组 (Parallel)":
            c_p1, c_p2 = st.columns(2)
            with c_p1:
                p_treatment_percent = st.number_input(
                    "试验组预期有效率 (%)", value=50.0, step=0.1,
                )
            with c_p2:
                p_control_percent = st.number_input(
                    "对照组预期有效率 (%)", value=40.0, step=0.1,
                )
            half_width_percent = st.number_input(
                "率差容许误差半宽 w (%)",
                value=10.0, step=0.5,
                help="两组率差 (pT − pC) 的置信区间半宽。",
            )
        else:
            p_treatment_percent = st.number_input(
                "预期有效率 P (%)", value=40.0, step=0.1,
            )
            p_control_percent = 0.0
            half_width_percent = st.number_input(
                "容许误差半宽 w (%)",
                value=10.0, step=0.5,
                help="95% CI 半宽不超过该值。",
            )
        return {
            "p_treatment_percent": p_treatment_percent,
            "p_control_percent": p_control_percent,
            "half_width_percent": half_width_percent,
            "margin_percent": 0.0, "direction": "higher",
        }

    # ---- 常规：独立两组 / 单臂 ----
    c_p1, c_p2 = st.columns(2)
    with c_p1:
        p_treatment_percent = st.number_input(
            "试验组预期有效率 (%)", value=60.0, step=0.1,
        )
    if group_type == "独立两组 (Parallel)":
        label_ctrl = "对照组预期有效率 (%)"
    else:
        label_ctrl = "历史/目标有效率 P0 (%)"
    with c_p2:
        p_control_percent = st.number_input(label_ctrl, value=40.0, step=0.1)

    # 方向选择（仅独立两组优效性时显示）
    if group_type == "独立两组 (Parallel)" and test_type == "优效性检验":
        direction_label = st.selectbox(
            "指标方向（优效）",
            ["高值更优（有效率更高）", "低值更优（事件率更低）"],
            index=0,
        )
        direction = "lower" if direction_label.startswith("低值更优") else "higher"
    else:
        direction = "higher"

    # 界值
    if test_type == "差异性检验":
        margin_percent = 0.0
    else:
        margin_label = (
            "非劣效界值 Δ" if test_type == "非劣效性检验"
            else "等效/优效界值 Δ"
        )
        margin_percent = st.number_input(
            f"{margin_label} (百分点)", value=5.0, step=0.1,
        )

    return {
        "p_treatment_percent": p_treatment_percent,
        "p_control_percent": p_control_percent,
        "margin_percent": margin_percent,
        "direction": direction,
    }


# ==============================================================================
# 主渲染入口
# ==============================================================================
def render(col_input, col_result, group_type: str, test_type: str,
           simon_design_type: str | None) -> None:
    variable_type = "分类变量 (率)"

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
        form_key = f"proportion_form_{group_type}_{test_type}"
        with st.form(key=form_key):
            st.subheader("临床预期", divider="gray")
            clinical = render_clinical_inputs(group_type, test_type,
                                              simon_design_type)

            # 解包 clinical 参数
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

            # 统计学设计参数
            stats = ui_common.render_statistics_design(
                variable_type, group_type, test_type,
            )
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

            # ---- 提交按钮 ----
            st.markdown("<br>", unsafe_allow_html=True)
            submitted = st.form_submit_button(
                label="GO: 开始计算样本量",
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

        elif (group_type == "独立两组 (Parallel)"
              and test_type == "精度分析"):
            pm = ("newcombe_wilson"
                  if precision_method.startswith("Newcombe") else "wald")
            prec = proportion_calculator.compute_two_arm_precision(
                p_t=p_treatment_percent / 100,
                p_c=p_control_percent / 100,
                r=alloc_ratio_val,
                ci_level=1 - alpha,
                half_width=half_width_percent / 100,
                method=pm,
            )
            nt, nc = int(prec["n_treatment"]), int(prec["n_control"])
            ntd = int(math.ceil(nt / (1 - dropout_percent / 100)))
            ncd = int(math.ceil(nc / (1 - dropout_percent / 100)))
            result = {
                "n_treatment": nt, "n_control": nc,
                "total_sample_size": nt + nc,
                "n_treatment_with_dropout": ntd,
                "n_control_with_dropout": ncd,
                "total_sample_size_with_dropout": ntd + ncd,
                "method": prec["method"], "precision": prec,
            }

        elif (group_type == "单臂设计 (Single Arm)"
              and test_type == "精度分析"):
            prec_method = (
                "clopper_pearson"
                if precision_method.startswith("确切") else "wald"
            )
            prec = proportion_calculator.compute_single_arm_precision(
                p=p_treatment_percent / 100,
                ci_level=1 - alpha,
                half_width=half_width_percent / 100,
                method=prec_method,
            )
            n = int(prec["n"])
            n_drop = int(math.ceil(n / (1 - dropout_percent / 100)))
            result = {
                "n_treatment": n, "n_control": 0,
                "total_sample_size": n,
                "n_treatment_with_dropout": n_drop,
                "n_control_with_dropout": 0,
                "total_sample_size_with_dropout": n_drop,
                "method": prec["method"], "precision": prec,
            }

        elif group_type == "配对两组 (Paired)":
            result = proportion_calculator.compute_paired_proportion(
                p_discordant_treat=p10_pct / 100.0,
                p_discordant_ctrl=p01_pct / 100.0,
                alpha=alpha,
                power=power_percent / 100.0,
                dropout_rate=dropout_percent / 100.0,
                test_type=test_type,
                margin=(margin_percent / 100.0
                        if test_type != "差异性检验" else 0.0),
                method=paired_rate_method,
            )

        elif group_type == "三组设计 (Three Arm)":
            if test_type == "差异性检验（率齐性）":
                result = proportion_calculator.compute_three_arm_proportion(
                    p1_pct / 100.0, p2_pct / 100.0, p3_pct / 100.0,
                    alpha=alpha, power=power_percent / 100.0,
                    dropout_rate=dropout_percent / 100.0,
                    method=three_arm_diff_method,
                )
            else:
                result = (
                    proportion_calculator.compute_three_arm_cochran_armitage(
                        p1_pct / 100.0, p2_pct / 100.0, p3_pct / 100.0,
                        alpha=alpha, power=power_percent / 100.0,
                        dropout_rate=dropout_percent / 100.0,
                        scores=(ca_t1, ca_t2, ca_t3),
                    )
                )

        else:
            # 独立两组 或 单臂（非精度分析）
            params = ProportionInput(
                p_treatment=p_treatment_percent / 100,
                p_control=p_control_percent / 100,
                alpha=alpha,
                power=power_percent / 100,
                test_type=test_type,
                allocation_ratio=alloc_ratio_val,
                dropout_rate=dropout_percent / 100,
                margin=margin_percent / 100,
                direction=(
                    direction
                    if group_type == "独立两组 (Parallel)"
                    else "higher"
                ),
                two_arm_diff_method=two_arm_diff_method,
                two_arm_margin_method=two_arm_margin_method,
            )
            if group_type == "独立两组 (Parallel)":
                result = proportion_calculator.compute_two_arm(params)
            else:
                result = proportion_calculator.compute_single_arm(
                    params, method=single_arm_method,
                )

    except Exception as e:
        error_message = str(e)

    # ------------------------------------------------------------------
    # 右侧：结果输出
    # ------------------------------------------------------------------
    with col_result:
        if error_message:
            st.error(f"计算冲突: {error_message}")
            return

        # ---------- Simon 二阶段 ----------
        if group_type == "Simon 二阶段" and simon_res:
            _render_simon_result(simon_res, p_control_percent,
                                 p_treatment_percent, alpha,
                                 power_percent, dropout_percent)
            return

        if not result:
            st.warning("无计算结果。")
            return

        total_n = result["total_sample_size_with_dropout"]
        n_treat = result.get("n_treatment_with_dropout", 0)
        n_ctrl = result.get("n_control_with_dropout", 0)
        theo_total = result["total_sample_size"]

        # ---------- 三组设计 ----------
        if group_type == "三组设计 (Three Arm)":
            _render_three_arm_result(result, test_type,
                                     p1_pct, p2_pct, p3_pct,
                                     alpha, power_percent,
                                     dropout_percent, alpha_side_label)
            return

        # ---------- 单臂 / 配对 ----------
        if group_type in ("单臂设计 (Single Arm)", "配对两组 (Paired)"):
            _render_single_or_paired_result(
                result, group_type, test_type,
                p_treatment_percent, p_control_percent,
                p10_pct, p01_pct, margin_percent,
                half_width_percent,
                alpha, power_percent, dropout_percent,
                alpha_side_label,
            )
            return

        # ---------- 独立两组 ----------
        _render_two_arm_result(
            result, test_type,
            p_treatment_percent, p_control_percent, margin_percent,
            alpha, power_percent, dropout_percent,
            alloc_ratio_str, alloc_ratio_val, alpha_side_label,
            total_n, n_treat, n_ctrl, theo_total,
        )


# ==============================================================================
# 各场景的输出渲染
# ==============================================================================

def _render_simon_result(simon_res, p0_pct, p1_pct,
                         alpha, power_pct, dropout_pct):
    """Simon 二阶段设计结果（Morandi KPI 卡片 + 极简版叙述）。"""
    design_label = "最优（optimal）" if simon_res["design_type"] == "optimal" else "极小极大（minimax）"

    # 1. Morandi KPI 卡片
    st.markdown(
        f"""
<div class="result-kpi-card">
    <div style="font-size:16px;font-weight:500;text-align:center;opacity:0.95;letter-spacing:0.5px;">
        Simon {simon_res['design_type'].capitalize()} 目标招募样本量
    </div>
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
        总计招募：{simon_res['total_sample_size_with_dropout']} 例
        （理论 {simon_res['n']} 例，已计入 {dropout_pct}% 脱落率）
    </div>
</div>
""",
        unsafe_allow_html=True,
    )

    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown('<div class="section-title">样本量计算方法</div>', unsafe_allow_html=True)

    # 2. 极简版叙述
    cn_text = (
        f"本研究主要终点为二分类结局（有效率），采用 Simon 两阶段{design_label}设计。"
        f"设定药物无效界值 P₀ = {p0_pct}%，预期目标有效率 P₁ = {p1_pct}%。"
        f"采用单侧检验，显著性水平 α = {alpha}，检验效能 1-β = {power_pct:.0f}%。"
        f"第一阶段入组 {simon_res['n1']} 例，若有效例数 ≤ {simon_res['r1']} 则提前终止"
        f"（PET = {simon_res['pet0'] * 100:.1f}%，在 P₀ 下）；"
        f"否则进入第二阶段，继续入组 {simon_res['n2']} 例。"
        f"最终若总有效例数 > {simon_res['r']}，则拒绝原假设。"
        f"考虑约 {dropout_pct:.0f}% 脱落率，实际计划最大入组 {simon_res['total_sample_size_with_dropout']} 例。"
    )

    st.markdown(
        f"""
<div class="statement-box">
    <p>{cn_text}</p>
</div>
""",
        unsafe_allow_html=True,
    )


def _render_three_arm_result(result, test_type,
                             p1_pct, p2_pct, p3_pct,
                             alpha, power_pct, dropout_pct,
                             alpha_side_label):
    """三组设计结果（Morandi KPI 卡片 + 极简版叙述）。"""
    theo_pg = int(result["n_per_group"])
    npg_drop = int(result["n_per_group_with_dropout"])
    theo_total = result["total_sample_size"]
    total_n = result["total_sample_size_with_dropout"]

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
        （每组 {npg_drop} 例，已计入 {dropout_pct}% 脱落率）
    </div>
</div>
""",
        unsafe_allow_html=True,
    )

    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown('<div class="section-title">样本量计算方法</div>', unsafe_allow_html=True)

    # 2. 极简版叙述
    if test_type == "差异性检验（率齐性）":
        cn_text = (
            f"本研究主要终点为二分类结局，设计为三组独立样本，"
            f"主要分析为率齐性检验（H₀: π₁ = π₂ = π₃）。"
            f"假定三组有效率分别约为 π₁ = {p1_pct:.1f}%、"
            f"π₂ = {p2_pct:.1f}%、π₃ = {p3_pct:.1f}%。"
            f"采用{alpha_side_label}，显著性水平 α = {alpha}，"
            f"检验效能 1-β = {power_pct:.0f}%，三组等额分配。"
            f"在上述参数条件下，依据 Pearson χ² 齐性检验（非中心 χ² 法）估算，"
            f"理论所需总样本量为 {theo_total} 例（每组 {theo_pg} 例）。"
            f"考虑约 {dropout_pct:.0f}% 脱落率，计划总入组 {total_n} 例。"
        )
    else:
        sc = result.get("ca_scores") or (0.0, 1.0, 2.0)
        cn_text = (
            f"本研究主要终点为二分类结局，设计为三组有序样本，"
            f"主要分析为 Cochran–Armitage 线性趋势检验（双侧）。"
            f"假定三组有效率分别约为 π₁ = {p1_pct:.1f}%、"
            f"π₂ = {p2_pct:.1f}%、π₃ = {p3_pct:.1f}%，"
            f"对应趋势分数为 ({sc[0]:g}, {sc[1]:g}, {sc[2]:g})。"
            f"采用{alpha_side_label}，显著性水平 α = {alpha}，"
            f"检验效能 1-β = {power_pct:.0f}%，三组等额分配。"
            f"在上述参数条件下，依据 Cochran–Armitage 趋势检验（非中心 χ² 法）估算，"
            f"理论所需总样本量为 {theo_total} 例（每组 {theo_pg} 例）。"
            f"考虑约 {dropout_pct:.0f}% 脱落率，计划总入组 {total_n} 例。"
        )

    st.markdown(
        f"""
<div class="statement-box">
    <p>{cn_text}</p>
</div>
""",
        unsafe_allow_html=True,
    )


def _render_single_or_paired_result(
    result, group_type, test_type,
    p_t_pct, p_c_pct, p10_pct, p01_pct, margin_pct,
    half_width_pct,
    alpha, power_pct, dropout_pct, alpha_side_label,
):
    """单臂或配对设计结果（Morandi KPI 卡片 + 极简版叙述）。"""
    total_n = result["total_sample_size_with_dropout"]
    theo_total = result["total_sample_size"]

    # 1. Morandi KPI 卡片（单臂/配对均用单值大卡片）
    kpi_title = "目标配对受试者数" if group_type == "配对两组 (Paired)" else "目标招募样本量"
    st.markdown(
        f"""
<div class="result-kpi-card">
    <div style="font-size:16px;font-weight:500;text-align:center;opacity:0.95;letter-spacing:0.5px;">
        {kpi_title}
    </div>
    <div style="text-align:center;font-size:48px;font-weight:700;margin:24px 0;">
        {total_n} 例
    </div>
    <div class="total-box">
        理论所需：{theo_total} 例（已计入 {dropout_pct}% 脱落率）
    </div>
</div>
""",
        unsafe_allow_html=True,
    )

    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown('<div class="section-title">样本量计算方法</div>', unsafe_allow_html=True)

    # 2. 极简版叙述
    if group_type == "配对两组 (Paired)":
        margin_txt = f"、界值 Δ = {margin_pct:.1f}%" if test_type != "差异性检验" else ""
        cn_text = (
            f"本研究主要终点为二分类结局（配对设计）。"
            f"假定不一致配对比例 π₁₀ ≈ {p10_pct:.2f}%、π₀₁ ≈ {p01_pct:.2f}%"
            f"（ψ = π₁₀ − π₀₁）{margin_txt}。"
            f"采用{alpha_side_label}，显著性水平 α = {alpha}，"
            f"检验效能 1-β = {power_pct:.0f}%。"
            f"在上述参数条件下，依据 {result.get('method', 'Nam Score 法')} 计算，"
            f"理论所需总样本量为 {theo_total} 例（对子）。"
            f"考虑约 {dropout_pct:.0f}% 脱落率，实际计划总入组 {total_n} 例。"
        )
    elif test_type == "精度分析":
        cn_text = (
            f"本研究为单臂率精度分析。"
            f"假定目标率约为 {p_t_pct:.1f}%，"
            f"设定 {(1 - alpha) * 100:.0f}% 置信区间半宽不超过 {half_width_pct:.1f}%。"
            f"在上述参数条件下，依据 {result.get('method', '正态近似法')} 计算，"
            f"理论所需样本量为 {theo_total} 例。"
            f"考虑约 {dropout_pct:.0f}% 脱落率，实际计划入组 {total_n} 例。"
        )
    else:
        margin_txt = f"、界值 Δ = {margin_pct:.1f}%" if test_type != "差异性检验" else ""
        cn_text = (
            f"本研究主要终点为二分类结局（单臂设计）。"
            f"根据临床预期，假定试验组有效率为 {p_t_pct:.1f}%，"
            f"历史对照目标值为 {p_c_pct:.1f}%{margin_txt}。"
            f"采用{alpha_side_label}，显著性水平 α = {alpha}，"
            f"检验效能 1-β = {power_pct:.0f}%。"
            f"在上述参数条件下，依据 {result.get('method', '正态近似法')} 计算，"
            f"理论所需总样本量为 {theo_total} 例。"
            f"考虑约 {dropout_pct:.0f}% 脱落率，实际计划总入组 {total_n} 例。"
        )

    st.markdown(
        f"""
<div class="statement-box">
    <p>{cn_text}</p>
</div>
""",
        unsafe_allow_html=True,
    )


def _render_two_arm_result(
    result, test_type,
    p_t_pct, p_c_pct, margin_pct,
    alpha, power_pct, dropout_pct,
    alloc_ratio_str, alloc_ratio_val, alpha_side_label,
    total_n, n_treat, n_ctrl, theo_total,
):
    """独立两组设计结果（保持原先精美的 Morandi 卡片设计，但极大简化文本，去除冗余描述）。"""

    # 1. 恢复原先精美的 Morandi KPI 容器大卡片
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
        总计招募：{total_n} 例（已计入 {dropout_pct}% 脱落率，分配比 {alloc_ratio_str}）
    </div>
</div>
""",
        unsafe_allow_html=True,
    )

    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown('<div class="section-title">样本量计算方法</div>', unsafe_allow_html=True)

    # 2. 极简版计算方法叙述（使用原版的 statement-box 样式，去除多余的 caption 和长段公式引用）
    margin_txt = f"、界值 Δ = {margin_pct:.1f}%" if test_type != "差异性检验" else ""
    cn_text = (
        f"本研究主要终点为二分类结局（有效率），采用独立两组{test_type[:2]}设计。"
        f"根据临床预期，试验组有效率预计为 {p_t_pct:.1f}%，对照组为 {p_c_pct:.1f}%{margin_txt}。"
        f"采用{alpha_side_label}，显著性水平 α = {alpha}，检验效能 1-β = {power_pct:.0f}%。"
        f"在上述参数条件下，依据 {result.get('method', '标准正态近似')} 法计算，理论所需总样本量为 {theo_total} 例。"
        f"考虑约 {dropout_pct:.0f}% 脱落率，实际计划最多入组 {total_n} 例（试验组 {n_treat} 例，对照组 {n_ctrl} 例）。"
    )

    st.markdown(
        f"""
<div class="statement-box">
    <p>{cn_text}</p>
</div>
""",
        unsafe_allow_html=True,
    )
