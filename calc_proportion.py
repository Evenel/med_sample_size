# calc_proportion.py
"""临床试验样本量计算 —— 分类变量（率）

覆盖设计：
  - 独立两组 (Parallel Two-Arm)：差异性 / 优效性 / 非劣效性 / 等效性 / 精度分析
  - 单臂设计 (Single Arm)：功效分析 / 精度分析
  - 配对两组 (Paired)：McNemar / Nam Score
  - 三组设计 (Three Arm)：率齐性 / Cochran–Armitage 趋势检验
  - Simon 二阶段 (Simon Two-Stage)

临床统计学参考：
  - Chow SC, Shao J, Wang H. Sample Size Calculations in Clinical Research (3rd Ed). CRC Press, 2017.
  - Farrington CP, Manning G. Test statistics and sample size formulae for comparative binomial
    trials with null hypothesis of non-zero risk difference or non-unity relative risk. Stat Med, 1990.
  - Newcombe RG. Interval estimation for the difference between independent proportions. Stat Med, 1998.
  - Nam J. Power and sample size requirements for non-inferiority in studies comparing two
    correlated proportions. Stat Med, 1997.
  - Simon R. Optimal two-stage designs for phase II clinical trials. Control Clin Trials, 1989.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Literal, Optional, Tuple

import numpy as np
from scipy.stats import beta, binom, binomtest, chi2, fisher_exact, ncx2, norm

# ---------------------------------------------------------------------------
# 类型别名
# ---------------------------------------------------------------------------
TestType = Literal["差异性检验", "优效性检验", "非劣效性检验", "等效性检验"]
Direction = Literal["higher", "lower"]
DiffMethod = Literal["auto", "pooled_z", "unpooled_z", "yates_correction", "fisher_exact"]
MarginMethod = Literal["auto", "fm_score", "wald"]

# ---------------------------------------------------------------------------
# 输入数据类
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ProportionInput:
    p_treatment: float
    p_control: Optional[float] = None
    alpha: float = 0.05
    power: float = 0.80
    test_type: TestType = "优效性检验"
    allocation_ratio: float = 1.0
    dropout_rate: float = 0.20
    margin: float = 0.0
    direction: Direction = "higher"
    two_arm_diff_method: Optional[str] = None   # auto / pooled_z / unpooled_z / yates_correction / fisher_exact
    two_arm_margin_method: Optional[str] = None  # auto / fm_score / wald


# ---------------------------------------------------------------------------
# 主计算器
# ---------------------------------------------------------------------------
class ProportionCalculator:
    """率（分类变量）样本量计算器。

    FM Score 使用 Farrington–Manning restricted MLE 求解 H₀ 边界下的约束极大似
    然估计，与 Wald（直接边界代入）在非零界值时结果不同，与 SAS/East/PASS 等
    商业软件可交叉验证。
    """

    # ---- 工具函数 ----------------------------------------------------------
    @staticmethod
    def _clip(p: float, lo: float = 1e-6, hi: float = 1.0 - 1e-6) -> float:
        return min(max(p, lo), hi)

    @staticmethod
    def _rates_in_normal_zone(pt: float, pc: float) -> bool:
        return 0.20 <= pt <= 0.80 and 0.20 <= pc <= 0.80

    @staticmethod
    def _expected_counts_ok(pt: float, pc: float, r: float, n_t: float) -> bool:
        if n_t <= 0:
            return False
        n_t_i = int(math.ceil(n_t))
        n_c_i = int(math.ceil(n_t_i * r))
        cells = (n_t_i * pt, n_t_i * (1 - pt), n_c_i * pc, n_c_i * (1 - pc))
        return min(cells) >= 5.0 - 1e-9

    def _z_alpha_beta(
        self, alpha: float, power: float, test: str
    ) -> Tuple[float, float]:
        """返回 (z_alpha, z_beta)，正确处理双侧/单侧/TOST。"""
        if test == "差异性检验":
            return float(norm.ppf(1 - alpha / 2.0)), float(norm.ppf(power))
        elif test == "等效性检验":
            # TOST: 每次单侧检验使用 α，power 在两侧各用 (1-β)/2 校正
            return float(norm.ppf(1 - alpha)), float(norm.ppf(1 - (1 - power) / 2.0))
        else:
            # 优效 / 非劣：单侧检验
            return float(norm.ppf(1 - alpha)), float(norm.ppf(power))

    # ---- Farrington–Manning Restricted MLE (核心) --------------------------
    def _fm_restricted_mle(
        self, pt: float, pc: float, delta: float, r: float
    ) -> Tuple[float, float, float]:
        """求解 H₀ 边界 p̃_T − p̃_C = δ 下的约束极大似然估计。

        三次方程在 (0,1) 内可能有两个根，选择使约束似然函数更大的根
        （与 Farrington–Manning 原文一致）。

        Returns:
            (p̃_C, p̃_T, V0) — 约束 MLE 下的对照组率、试验组率、合并零假设方差。
        """
        # Farrington–Manning 三次方程系数（以 p̃_C 为未知数）
        a = 1.0 + r
        b = -(1.0 + r + pt + r * pc + delta * (r + 2.0))
        c = delta ** 2 + delta * (2.0 * pt + r + 1.0) + pt + r * pc
        d = -pt * delta * (1.0 + delta)

        # 数值求解三次方程：取 (0,1) 内的所有实根
        roots = np.roots([a, b, c, d])
        candidates: List[float] = []
        for root in roots:
            if np.isreal(root):
                rv = float(root.real)
                p_t = rv + delta
                if 1e-8 < rv < 1.0 - 1e-8 and 1e-8 < p_t < 1.0 - 1e-8:
                    candidates.append(rv)

        if not candidates:
            # 退化：边界直接代入（Wald 型后备）
            x = (pt + pc - delta) / 2.0
            p_tilde_c = self._clip(x, 0.001, 0.999)
        elif len(candidates) == 1:
            p_tilde_c = candidates[0]
        else:
            # 多个候选：选择约束对数似然较大的根
            def constrained_loglik(pc_: float) -> float:
                pt_ = pc_ + delta
                ll = 0.0
                if 0 < pt_ < 1:
                    ll += pt * math.log(pt_) + (1 - pt) * math.log(1 - pt_)
                else:
                    return -float("inf")
                if 0 < pc_ < 1:
                    ll += r * pc * math.log(pc_) + r * (1 - pc) * math.log(1 - pc_)
                else:
                    return -float("inf")
                return ll

            p_tilde_c = max(candidates, key=constrained_loglik)

        p_tilde_t = self._clip(p_tilde_c + delta, 0.001, 0.999)
        p_tilde_c = self._clip(p_tilde_c, 0.001, 0.999)

        V0 = p_tilde_t * (1 - p_tilde_t) + p_tilde_c * (1 - p_tilde_c) / r
        return p_tilde_c, p_tilde_t, float(V0)

    def _wald_boundary_v0(
        self, pt: float, pc: float, delta: float, r: float
    ) -> float:
        """Wald 型 H₀ 边界方差：直接在边界 p_T = p_C + δ 处计算（非约束 MLE）。

        FM Score 与 Wald 的核心区别：
        - FM: 使用 restricted MLE（三次方程解）得到 p̃_C, p̃_T
        - Wald: 直接令 p_T0 = min(max(pc + δ, ε), 1−ε)，p_C0 = pc
        """
        pt0 = self._clip(pc + delta, 0.001, 0.999)
        return float(pt0 * (1 - pt0) + pc * (1 - pc) / r)

    # ---- 两独立样本率差样本量 ------------------------------------------------
    def compute_two_arm(self, params: ProportionInput) -> Dict[str, Any]:
        pt, pc, r = params.p_treatment, params.p_control, params.allocation_ratio
        margin, test = params.margin, params.test_type
        direction = params.direction

        # ---------- 差异性检验 ----------
        if test == "差异性检验":
            return self._compute_two_arm_difference(
                pt, pc, r, params.alpha, params.power,
                params.dropout_rate, params.two_arm_diff_method,
            )

        # ---------- 优效 / 非劣 / 等效 ----------
        return self._compute_two_arm_margin(
            pt, pc, r, margin, test, direction,
            params.alpha, params.power, params.dropout_rate,
            params.two_arm_margin_method,
        )

    # ---- 差异性检验：多种方法 ------------------------------------------------
    def _compute_two_arm_difference(
        self,
        pt: float, pc: float, r: float,
        alpha: float, power: float, dropout: float,
        method: Optional[str],
    ) -> Dict[str, Any]:
        """独立性两组率差异性检验（双侧，H₀: pT = pC）。

        提供方法：
          auto           — 自动选择（期望频数≥5 → pooled_z，否则 → fisher_exact）
          pooled_z       — 合并方差正态近似（与 Pearson χ² 等价，药政常用基准）
          unpooled_z     — 非合并方差简化公式（国内方案常用，公式简洁）
          yates_correction — 含 Yates 连续性校正（小样本保守选择）
          fisher_exact   — Fisher 确切检验（枚举功效，控制 I 类错误）
        """
        z_alpha, z_beta = self._z_alpha_beta(alpha, power, "差异性检验")
        effect = abs(pt - pc)
        if effect < 1e-12:
            raise ValueError("两组预期有效率相同，效应量为零，无法计算样本量。")

        V1 = pt * (1 - pt) + pc * (1 - pc) / r
        normal_zone = self._rates_in_normal_zone(pt, pc)

        # ---- 方法选择 ----
        dm = method or "auto"
        if dm == "auto":
            p_avg = (pt + r * pc) / (1 + r)
            V0_pooled = p_avg * (1 - p_avg) * (1 + 1 / r)
            nt_plan = ((z_alpha * math.sqrt(V0_pooled) + z_beta * math.sqrt(V1)) ** 2) / (effect ** 2)
            dm = "pooled_z" if self._expected_counts_ok(pt, pc, r, nt_plan) else "fisher_exact"

        # ---- Fisher 精确检验 ----
        if dm == "fisher_exact":
            return self._fisher_exact_sample_size_diff(pt, pc, r, alpha, power, dropout)

        # ---- 非合并方差 Z ----
        if dm == "unpooled_z":
            # 非合并方差：H₀ 与 H₁ 使用相同方差 V1
            nt = ((z_alpha + z_beta) ** 2 * V1) / (effect ** 2)
            out = self._format_result(nt, nt * r, dropout)
            out["method"] = "两独立样本率差：非合并方差正态近似"
            out["method_ref"] = (
                "Chow SC, et al. Sample Size Calculations in Clinical Research (3rd Ed), "
                "Section 4.1.1 (Two-Sample Z-Test with Unpooled Variance)."
            )
            out["rate_zone"] = "normal" if normal_zone else "extreme"
            out["rate_note"] = (
                "非合并方差公式（两样本 Z 检验简化形式），H₀ 与 H₁ 使用相同方差估计；"
                "国内临床试验方案常用。与合并方差法相比，在等比例分配时样本量略小。"
                + (
                    " 当前至少一组率超出 [0.2,0.8]，正态近似可能偏倚，"
                    "建议分析阶段使用 Fisher 精确检验或 CMH 检验进行敏感性分析。"
                    if not normal_zone
                    else ""
                )
            )
            return out

        # ---- Yates 连续性校正 ----
        if dm == "yates_correction":
            p_avg = (pt + r * pc) / (1 + r)
            V0 = p_avg * (1 - p_avg) * (1 + 1 / r)
            # Fleiss 连续性校正公式 (Fleiss, 1981)
            n_prime = 1 + r
            nt_raw = (
                n_prime
                / (4 * r * effect**2)
                * (
                    1.0
                    + math.sqrt(
                        1.0
                        + 4 * r * effect
                        / (n_prime * (z_alpha * math.sqrt(V0) + z_beta * math.sqrt(V1)) ** 2)
                    )
                )
                * (z_alpha * math.sqrt(V0) + z_beta * math.sqrt(V1)) ** 2
            )
            out = self._format_result(nt_raw, nt_raw * r, dropout)
            out["method"] = "两独立样本率差：Pearson χ² 含连续性校正（Yates 校正）"
            out["method_ref"] = (
                "Fleiss JL. Statistical Methods for Rates and Proportions (2nd Ed). Wiley, 1981, "
                "Section 3.3."
            )
            out["rate_zone"] = "yates"
            out["rate_note"] = (
                "Yates 连续性校正使 I 类错误更保守；小样本（期望频数 <5）时推荐，"
                "但检验效能略低于无校正版本。注册试验可在方案中说明使用该保守估计。"
            )
            return out

        # ---- 合并方差 Z（pooled_z，默认推荐） ----
        p_avg = (pt + r * pc) / (1 + r)
        V0 = p_avg * (1 - p_avg) * (1 + 1 / r)
        nt_raw = ((z_alpha * math.sqrt(V0) + z_beta * math.sqrt(V1)) ** 2) / (effect ** 2)
        nt_int = int(math.ceil(nt_raw))
        nc_int = int(math.ceil(nt_int * r))
        exp_ok = self._expected_counts_ok(pt, pc, r, nt_raw)

        if exp_ok and normal_zone:
            note = (
                "合并方差正态近似（与 Pearson χ² 独立性检验等价）。"
                "按当前假设估计样本量下，四格表期望频数均 ≥5（Cochran 规则），正态近似可靠。"
            )
            zone = "normal"
        elif exp_ok:
            note = (
                "合并方差正态近似；至少一组率不在 [0.2,0.8]，"
                "但按当前假设估计的样本量下期望频数仍 ≥5，正态近似仍可接受。"
                "确证性分析建议以 CMH 检验或 Fisher 精确检验作为敏感性分析。"
            )
            zone = "extreme_ok"
        else:
            note = (
                "合并方差正态近似；按当前假设估计样本量下部分期望频数 <5，"
                "正态近似可能偏倚，分析阶段建议使用 Fisher 精确检验或改用连续性校正法。"
            )
            zone = "extreme_sparse"

        out = self._format_result(nt_raw, nt_raw * r, dropout)
        out["method"] = "两独立样本率差：合并方差正态近似（与 Pearson χ² 独立性检验等价）"
        out["method_ref"] = (
            "Chow SC, et al. Sample Size Calculations in Clinical Research (3rd Ed), "
            "Section 4.1.1 (Two-Sample Z-Test with Pooled Variance)."
        )
        out["rate_zone"] = zone
        out["rate_note"] = note
        out["expected_min_cell"] = float(
            min(nt_int * pt, nt_int * (1 - pt), nc_int * pc, nc_int * (1 - pc))
        )
        return out

    # ---- 优效 / 非劣 / 等效：FM Score / Wald --------------------------------
    def _compute_two_arm_margin(
        self,
        pt: float, pc: float, r: float,
        margin: float, test: str, direction: str,
        alpha: float, power: float, dropout: float,
        method: Optional[str],
    ) -> Dict[str, Any]:
        """独立两组率优效/非劣/等效检验样本量。

        提供方法：
          auto / fm_score — Farrington–Manning Score（restricted MLE，药政金标准）
          wald            — Wald 型边界直接代入（大样本简化）
        """
        z_alpha, z_beta = self._z_alpha_beta(alpha, power, test)
        V1 = pt * (1 - pt) + pc * (1 - pc) / r

        # ---- 确定效应量 δ₀ 与检验界值 ----
        if test == "优效性检验":
            if direction == "higher":
                effect = (pt - pc) - margin  # H₁: pT − pC > Δ
                boundary_delta = margin
            else:  # lower: 事件率更低更优
                effect = (pc - pt) - margin  # H₁: pC − pT > Δ
                boundary_delta = -margin
        elif test == "非劣效性检验":
            effect = (pt - pc) + margin  # H₁: pT − pC > −Δ
            boundary_delta = -margin
        elif test == "等效性检验":
            # TOST: 取两个单侧中较难满足的一侧（效应量较小者）
            effect_upper = margin - (pt - pc)
            effect_lower = margin + (pt - pc)
            effect = min(effect_upper, effect_lower)
            boundary_delta = margin  # 等效两侧界值对称
        else:
            raise ValueError(f"不支持的检验类型: {test}")

        if effect <= 0:
            raise ValueError(
                f"当前设定下效应量 ≤ 0，无法满足【{test}】要求。"
                f"请调整预期率、界值 Δ 或检验方向。"
            )

        normal_zone = self._rates_in_normal_zone(pt, pc)

        # ---- 方法选择 ----
        mm = method or "auto"
        if mm == "auto":
            mm = "fm_score"

        # ---- Wald 法 ----
        if mm == "wald":
            if test == "等效性检验":
                # TOST: 取两侧 Wald 方差较大者
                delta_sup = boundary_delta
                delta_inf = -boundary_delta
                V0_wald = max(
                    self._wald_boundary_v0(pt, pc, delta_sup, r),
                    self._wald_boundary_v0(pt, pc, delta_inf, r),
                )
            else:
                V0_wald = self._wald_boundary_v0(pt, pc, boundary_delta, r)

            nt = ((z_alpha * math.sqrt(V0_wald) + z_beta * math.sqrt(V1)) ** 2) / (effect ** 2)
            out = self._format_result(nt, nt * r, dropout)

            if test == "等效性检验":
                out["method"] = "两独立样本率：Wald TOST（双单侧正态近似）"
            else:
                out["method"] = "两独立样本率差：Wald 正态近似（H₀ 边界直接代入）"
            out["method_ref"] = (
                "Chow SC, et al. Sample Size Calculations in Clinical Research (3rd Ed), "
                "Section 4.2 (Two-Sample Z-Test with Unequal Variance)."
            )
            out["rate_zone"] = "normal" if normal_zone else "extreme"
            out["rate_note"] = (
                "Wald 法直接以 H₀ 边界值代入方差计算，不进行约束 MLE 迭代。"
                "大样本下与 FM Score 接近，但界值较大或率极端时可能偏保守。"
                "注册试验推荐以 FM Score 为主，Wald 为敏感性分析。"
            )
            return out

        # ---- FM Score（默认推荐） ----
        if test == "等效性检验":
            # TOST: 分别求两侧边界 restricted MLE 方差，取大者（保守）
            _, _, V0_sup = self._fm_restricted_mle(pt, pc, boundary_delta, r)
            _, _, V0_inf = self._fm_restricted_mle(pt, pc, -boundary_delta, r)
            V0_fm = max(V0_sup, V0_inf)
        else:
            _, _, V0_fm = self._fm_restricted_mle(pt, pc, boundary_delta, r)

        nt = ((z_alpha * math.sqrt(V0_fm) + z_beta * math.sqrt(V1)) ** 2) / (effect ** 2)
        out = self._format_result(nt, nt * r, dropout)

        if test == "等效性检验":
            out["method"] = "两独立样本率：Farrington–Manning TOST（双单侧 FM Score）"
            out["rate_note"] = (
                "等效性检验采用双单侧检验（TOST）框架，每次单侧检验使用显著性水平 α。"
                "零假设 H₀₁: pT − pC ≤ −Δ，H₀₂: pT − pC ≥ +Δ。"
                "FM Score 使用约束极大似然估计（restricted MLE）求解 H₀ 边界方差，"
                "与 SAS PROC POWER、PASS、East 等商业软件方法一致。"
            )
        else:
            test_label = {"优效性检验": "优效性", "非劣效性检验": "非劣效性"}.get(test, test[:2])
            direction_label = {
                "优效性检验": "优效" if direction == "higher" else "优效（低值更优）",
                "非劣效性检验": "非劣效",
            }.get(test, test[:2])
            out["method"] = f"两独立样本率差：Farrington–Manning Score 检验（{test_label}）"
            out["rate_note"] = (
                f"采用{test_label}检验（单侧 α={alpha}）。"
                "Farrington–Manning Score 使用 H₀ 边界下的约束极大似然估计"
                "（restricted MLE）计算检验统计量的方差，"
                "为国际多中心临床试验及药政申报（FDA/EMA/NMPA）推荐的基准方法之一。"
                + (
                    " 当前至少一组率在 [0.2,0.8] 之外，确证性试验建议以无条件精确检验"
                    "（如 Chan 法，需 StatXact 或 Cytel East 实现）进行复核。"
                    if not normal_zone
                    else ""
                )
            )
        out["method_ref"] = (
            "Farrington CP, Manning G. Test statistics and sample size formulae for "
            "comparative binomial trials with null hypothesis of non-zero risk difference "
            "or non-unity relative risk. Stat Med, 1990; 9: 1447–1454."
        )
        out["rate_zone"] = "normal" if normal_zone else "extreme"
        return out

    # ---- Fisher 精确检验（差异性）--------------------------------------------
    @staticmethod
    def _fisher_two_sample_power(
        n1: int, n2: int, p1: float, p2: float, alpha: float
    ) -> float:
        """双侧 Fisher 精确检验在 H₀: p1 = p2 下的近似功效（枚举所有 2×2 表）。"""
        s = 0.0
        for a in range(n1 + 1):
            pmf_a = float(binom.pmf(a, n1, p1))
            if pmf_a < 1e-15:
                continue
            for c in range(n2 + 1):
                pmf_c = float(binom.pmf(c, n2, p2))
                if pmf_c < 1e-20:
                    continue
                b, d = n1 - a, n2 - c
                try:
                    pv = float(fisher_exact([[a, b], [c, d]], alternative="two-sided").pvalue)
                except ValueError:
                    continue
                if pv <= alpha:
                    s += pmf_a * pmf_c
        return float(s)

    def _fisher_exact_sample_size_diff(
        self, pt: float, pc: float, r: float,
        alpha: float, power: float, dropout: float,
    ) -> Dict[str, Any]:
        """差异性：Fisher 精确检验枚举搜索最小样本量。

        若正态近似估算样本量 > 50，自动降级为 Yates 连续性校正，
        以避免三层枚举导致 CPU 假死。
        """
        z_a, z_b = self._z_alpha_beta(alpha, power, "差异性检验")
        effect = abs(pt - pc)
        p_avg = (pt + r * pc) / (1 + r)
        V0 = p_avg * (1 - p_avg) * (1 + 1 / r)
        V1 = pt * (1 - pt) + pc * (1 - pc) / r
        n0 = int(max(8, math.ceil(((z_a * math.sqrt(V0) + z_b * math.sqrt(V1)) ** 2) / (effect ** 2))))

        # ---- 安全保护：正态近似估算样本量 > 50 时禁止 Fisher 枚举，降级为 Yates ----
        FISHER_SAFE_THRESHOLD = 50
        if n0 > FISHER_SAFE_THRESHOLD:
            # 降级为 Yates 连续性校正（Fleiss 公式）
            n_prime = 1 + r
            nt_raw = (
                n_prime
                / (4 * r * effect**2)
                * (
                    1.0
                    + math.sqrt(
                        1.0
                        + 4 * r * effect
                        / (n_prime * (z_a * math.sqrt(V0) + z_b * math.sqrt(V1)) ** 2)
                    )
                )
                * (z_a * math.sqrt(V0) + z_b * math.sqrt(V1)) ** 2
            )
            out = self._format_result(nt_raw, nt_raw * r, dropout)
            out["method"] = (
                "两独立样本率差：Pearson χ² 含连续性校正（Yates 校正）"
                "【由 Fisher 精确检验自动降级】"
            )
            out["method_ref"] = (
                "Fleiss JL. Statistical Methods for Rates and Proportions (2nd Ed). Wiley, 1981, "
                "Section 3.3."
            )
            out["rate_zone"] = "yates_fallback"
            out["rate_note"] = (
                f"由于正态近似估算样本量 n₀≈{n0} > {FISHER_SAFE_THRESHOLD}，"
                "为防止 Fisher 精确枚举计算假死，已自动切换为连续性校正 Z 检验"
                "（Yates 校正）。Yates 校正在小样本时较保守，检验效能略低于无校正版本；"
                "但在此样本量下与 Fisher 确切检验结果差异极小，可作为注册试验的保守估计。"
            )
            return out

        n_start = max(5, n0 - 30)
        n_max = 8000
        for n1 in range(n_start, n_max + 1):
            n2 = int(math.ceil(n1 * r))
            pw = self._fisher_two_sample_power(n1, n2, pt, pc, alpha)
            if pw >= power:
                out = self._format_result(float(n1), float(n2), dropout)
                out["method"] = "两独立样本率差：Fisher 确切概率法（双侧，枚举功效）"
                out["method_ref"] = (
                    "Fisher RA. The Design of Experiments (8th Ed). Oliver & Boyd, 1966."
                )
                out["rate_zone"] = "fisher"
                out["rate_note"] = (
                    "基于 Fisher 确切检验（超几何分布）枚举所有可能的 2×2 表计算功效，"
                    "不依赖正态近似假设。适用于小样本或期望频数 <5 的场景。"
                    "大样本时结果与 Pearson χ² 检验（合并方差）趋于一致。"
                )
                out["actual_power_at_alt"] = float(pw)
                return out
        raise ValueError(
            "Fisher 精确检验在可搜索范围内未找到满足效能的样本量；请放宽效应量或改用正态近似法。"
        )

    # ---- 单臂率样本量（功效分析）---------------------------------------------
    @staticmethod
    def _cp_ci(x: int, n: int, alpha: float) -> Tuple[float, float]:
        """Clopper–Pearson 精确置信区间。"""
        if n <= 0:
            raise ValueError("n 必须为正整数。")
        lo = 0.0 if x == 0 else float(beta.ppf(alpha / 2.0, x, n - x + 1))
        hi = 1.0 if x == n else float(beta.ppf(1 - alpha / 2.0, x + 1, n - x))
        return lo, hi

    @staticmethod
    def _wilson_ci(x: int, n: int, alpha: float) -> Tuple[float, float]:
        """Wilson Score 置信区间。"""
        if n <= 0:
            raise ValueError("n 必须为正整数。")
        z = float(norm.ppf(1 - alpha / 2.0))
        phat = x / n
        denom = 1.0 + z ** 2 / n
        center = (phat + z ** 2 / (2 * n)) / denom
        half = (z / denom) * math.sqrt((phat * (1 - phat) + z ** 2 / (4 * n)) / n)
        return max(0.0, center - half), min(1.0, center + half)

    def compute_single_arm(
        self, params: ProportionInput,
        method: Literal["score", "exact"] = "score",
    ) -> Dict[str, Any]:
        """单臂率样本量（功效分析）。

        H₀ 边界点 p0_bound：
          - 差异性：p0_bound = p0（双侧）
          - 优效性：p0_bound = p0 + Δ（单侧）
          - 非劣效性：p0_bound = p0 − Δ（单侧）
          - 等效性：TOST，双侧边界
        """
        pt, p0, margin, test = params.p_treatment, params.p_control, params.margin, params.test_type
        if p0 is None:
            raise ValueError("单臂设计需要提供历史/目标有效率 P0。")
        if not (0 < pt < 1 and 0 < p0 < 1):
            raise ValueError("P0/P1 必须位于 (0,1)。")

        if test == "差异性检验":
            alpha_side = "two-sided"
            p0_bound = p0
            z_alpha = float(norm.ppf(1 - params.alpha / 2.0))
            z_beta = float(norm.ppf(params.power))
            effect = abs(pt - p0)
        elif test == "优效性检验":
            alpha_side = "one-sided"
            p0_bound = p0 + margin
            z_alpha = float(norm.ppf(1 - params.alpha))
            z_beta = float(norm.ppf(params.power))
            effect = pt - (p0 + margin)
        elif test == "非劣效性检验":
            alpha_side = "one-sided"
            p0_bound = p0 - margin
            z_alpha = float(norm.ppf(1 - params.alpha))
            z_beta = float(norm.ppf(params.power))
            effect = pt - (p0 - margin)
        elif test == "等效性检验":
            alpha_side = "tost"
            p0_bound = p0
            z_alpha = float(norm.ppf(1 - params.alpha))
            z_beta = float(norm.ppf(1 - (1 - params.power) / 2.0))
            effect = margin - abs(pt - p0)
        else:
            raise ValueError(f"不支持的检验类型: {test}")

        if effect <= 0:
            raise ValueError(
                f"当前设定下效应量 ≤ 0，无法满足【{test}】要求。"
                f"请调整 P0/P1 或界值 Δ。"
            )

        p0_bound = self._clip(p0_bound, 0.001, 0.999)

        if method == "score":
            n = (
                (z_alpha * math.sqrt(p0_bound * (1 - p0_bound))
                 + z_beta * math.sqrt(pt * (1 - pt))) ** 2
            ) / (effect ** 2)
            out = self._format_result(n, 0, params.dropout_rate)

            test_label_map = {
                "差异性检验": "差异性（双侧）",
                "优效性检验": "优效性（单侧）",
                "非劣效性检验": "非劣效性（单侧）",
                "等效性检验": "等效性（TOST，单侧 α）",
            }
            out["method"] = f"单臂率：Score 正态近似（{test_label_map.get(test, test)}）"
            out["method_ref"] = (
                "Chow SC, et al. Sample Size Calculations in Clinical Research (3rd Ed), "
                "Section 4.4 (Single-Arm Trials)."
            )
            out["alpha_side"] = alpha_side
            out["p0_bound"] = float(p0_bound)
            return out

        # ---- 精确二项 ----
        if method != "exact":
            raise ValueError("method 仅支持 'score' 或 'exact'。")

        n0 = int(
            max(
                5,
                math.ceil(
                    ((z_alpha * math.sqrt(p0_bound * (1 - p0_bound))
                      + z_beta * math.sqrt(pt * (1 - pt))) ** 2) / (effect ** 2)
                ),
            )
        )
        n_start, n_max = max(5, n0 - 50), 5000

        # 缓存
        sf_cache: Dict[Tuple[int, float, int], float] = {}
        cdf_cache: Dict[Tuple[int, float, int], float] = {}

        def sf(n_: int, p_: float, k_: int) -> float:
            key = (n_, p_, k_)
            if key not in sf_cache:
                sf_cache[key] = float(binom.sf(k_ - 1, n_, p_))
            return sf_cache[key]

        def cdf(n_: int, p_: float, k_: int) -> float:
            key = (n_, p_, k_)
            if key not in cdf_cache:
                cdf_cache[key] = float(binom.cdf(k_, n_, p_))
            return cdf_cache[key]

        def exact_p_two_sided(n_: int, p_: float, x_: int) -> float:
            phat = x_ / n_
            return min(2.0 * (sf(n_, p_, x_) if phat >= p_ else cdf(n_, p_, x_)), 1.0)

        for n in range(n_start, n_max + 1):
            reject = [False] * (n + 1)
            if test == "差异性检验":
                for x in range(n + 1):
                    if exact_p_two_sided(n, p0_bound, x) <= params.alpha:
                        reject[x] = True
            elif test in ("优效性检验", "非劣效性检验"):
                for x in range(n + 1):
                    if sf(n, p0_bound, x) <= params.alpha:
                        reject[x] = True
            elif test == "等效性检验":
                p_lo = self._clip(p0 - margin, 0.001, 0.999)
                p_hi = self._clip(p0 + margin, 0.001, 0.999)
                for x in range(n + 1):
                    if sf(n, p_lo, x) <= params.alpha and cdf(n, p_hi, x) <= params.alpha:
                        reject[x] = True
            else:
                raise ValueError(f"不支持的检验类型: {test}")

            # Power & actual alpha
            pmf_h1 = [float(binom.pmf(x, n, pt)) for x in range(n + 1)]
            pmf_h0 = [float(binom.pmf(x, n, p0_bound)) for x in range(n + 1)]
            pw = sum(pmf_h1[x] for x in range(n + 1) if reject[x])
            actual_alpha = sum(pmf_h0[x] for x in range(n + 1) if reject[x])
            if pw >= params.power:
                out = self._format_result(n, 0, params.dropout_rate)
                out["method"] = "单臂率：精确二项检验（Exact Binomial，枚举功效）"
                out["method_ref"] = (
                    "Chow SC, et al. Sample Size Calculations in Clinical Research (3rd Ed), "
                    "Section 4.4.2."
                )
                out["alpha_side"] = alpha_side
                out["p0_bound"] = float(p0_bound)
                out["actual_power_at_p1"] = float(pw)
                out["actual_alpha_at_p0"] = float(actual_alpha)
                return out
        raise ValueError("精确二项在搜索上限内未找到满足效能的样本量，请检查参数或改用近似法。")

    # ---- 单臂率精度分析 ------------------------------------------------------
    def compute_single_arm_precision(
        self,
        p: float,
        ci_level: float,
        half_width: float,
        method: Literal["clopper_pearson", "wilson", "wald"] = "clopper_pearson",
        strategy: Literal["point"] = "point",
        n_max: int = 200000,
    ) -> Dict[str, Any]:
        """单臂率精度分析：给定目标置信区间半宽，反推最小样本量。"""
        if not (0 < p < 1):
            raise ValueError("p 必须位于 (0,1)。")
        alpha = 1.0 - ci_level

        if method == "wald":
            z = float(norm.ppf(1 - alpha / 2.0))
            n = max(1, int(math.ceil(z ** 2 * p * (1 - p) / (half_width ** 2))))
            hw = z * math.sqrt(p * (1 - p) / n)
            return {
                "n": n, "x_reference": int(round(n * p)),
                "p_reference": float(p), "ci_level": float(ci_level),
                "half_width_target": float(half_width),
                "half_width_achieved": float(hw),
                "ci_lower": float(p - hw), "ci_upper": float(p + hw),
                "method": "单臂率精度分析：正态近似法（Wald）",
            }

        z = float(norm.ppf(1 - alpha / 2.0))
        n0 = int(math.ceil(z ** 2 * p * (1 - p) / (half_width ** 2)))
        n = max(5, n0 - 50)

        def interval(n_: int) -> Tuple[float, float, int]:
            x = int(round(n_ * p))
            if method == "clopper_pearson":
                lo, hi = self._cp_ci(x, n_, alpha)
            elif method == "wilson":
                lo, hi = self._wilson_ci(x, n_, alpha)
            else:
                raise ValueError("method 仅支持 clopper_pearson / wilson / wald")
            return lo, hi, x

        while n <= n_max:
            lo, hi, x = interval(n)
            hw = (hi - lo) / 2.0
            if hw <= half_width:
                method_name = {
                    "clopper_pearson": "Clopper–Pearson 精确置信区间",
                    "wilson": "Wilson Score 置信区间",
                }.get(method, method)
                return {
                    "n": n, "x_reference": x, "p_reference": float(p),
                    "ci_level": float(ci_level),
                    "half_width_target": float(half_width),
                    "half_width_achieved": float(hw),
                    "ci_lower": float(lo), "ci_upper": float(hi),
                    "method": f"单臂率精度分析：{method_name}",
                }
            n += 1
        raise ValueError("在搜索上限内未找到满足精度要求的样本量，请检查参数。")

    # ---- 两臂率差精度分析 ----------------------------------------------------
    def compute_two_arm_precision(
        self,
        p_t: float, p_c: float, r: float,
        ci_level: float, half_width: float,
        method: Literal["wald", "newcombe_wilson"] = "wald",
        n_max: int = 250000,
    ) -> Dict[str, Any]:
        """两臂率差 (p_T − p_C) 置信区间半宽精度下的最小样本量。"""
        if not (0 < p_t < 1 and 0 < p_c < 1):
            raise ValueError("两组预期率须位于 (0,1)。")
        alpha = 1.0 - ci_level
        z = float(norm.ppf(1 - alpha / 2.0))

        if method == "wald":
            n_t = max(1, int(math.ceil(z ** 2 * (p_t * (1 - p_t) + p_c * (1 - p_c) / r) / (half_width ** 2))))
            n_c = int(math.ceil(n_t * r))
            se = math.sqrt(p_t * (1 - p_t) / n_t + p_c * (1 - p_c) / n_c)
            hw = z * se
            return {
                "n_treatment": n_t, "n_control": n_c,
                "total_sample_size": n_t + n_c,
                "ci_level": float(ci_level),
                "half_width_target": float(half_width),
                "half_width_achieved": float(hw),
                "method": "两臂率差精度：Wald 正态近似（率差标准误）",
                "p_treatment": float(p_t), "p_control": float(p_c),
            }

        # Newcombe-Wilson (method 10)
        for n_t in range(5, n_max + 1):
            n_c = int(math.ceil(n_t * r))
            x_t = min(max(int(round(n_t * p_t)), 0), n_t)
            x_c = min(max(int(round(n_c * p_c)), 0), n_c)
            lo_t, hi_t = self._wilson_ci(x_t, n_t, alpha)
            lo_c, hi_c = self._wilson_ci(x_c, n_c, alpha)
            diff = p_t - p_c
            inner = max(0.0, (p_t - lo_t) ** 2 + (hi_c - p_c) ** 2)
            outer = max(0.0, (hi_t - p_t) ** 2 + (p_c - lo_c) ** 2)
            L = diff - math.sqrt(inner)
            U = diff + math.sqrt(outer)
            asym = max(abs(U - diff), abs(diff - L))
            if asym <= half_width:
                return {
                    "n_treatment": n_t, "n_control": n_c,
                    "total_sample_size": n_t + n_c,
                    "ci_level": float(ci_level),
                    "half_width_target": float(half_width),
                    "half_width_achieved": float(asym),
                    "ci_lower_diff": float(L), "ci_upper_diff": float(U),
                    "method": "两臂率差精度：Newcombe–Wilson Score 法（方法 10）",
                    "p_treatment": float(p_t), "p_control": float(p_c),
                }
        raise ValueError("在搜索上限内未找到满足精度要求的样本量，请检查参数或改用 Wald。")

    # ---- Simon 二阶段设计 ----------------------------------------------------
    def compute_simon_two_stage(
        self,
        p0: float, p1: float,
        alpha: float, power: float, dropout_rate: float,
        design_type: Literal["optimal", "minimax"] = "optimal",
        n_max: int = 100,
    ) -> dict:
        """精确 Simon 两阶段单臂设计（枚举搜索）。

        - optimal: 在可行设计中最小化 EN(p0)
        - minimax: 在可行设计中最小化总样本量 n
        """
        if not (0 < p0 < p1 < 1):
            raise ValueError("须满足 0 < P0 < P1 < 1。")
        if not (0 < alpha < 1 and 0 < power < 1):
            raise ValueError("alpha/power 必须位于 (0,1)。")

        pmf_cache, tail_cache = {}, {}

        def get_pmf(n1: int, p: float) -> list:
            key = (n1, p)
            if key not in pmf_cache:
                pmf_cache[key] = [float(binom.pmf(x, n1, p)) for x in range(n1 + 1)]
            return pmf_cache[key]

        def get_tail(n2: int, p: float) -> list:
            key = (n2, p)
            if key not in tail_cache:
                t = [1.0] * (n2 + 2)
                for k in range(1, n2 + 1):
                    t[k] = float(binom.sf(k - 1, n2, p))
                t[n2 + 1] = 0.0
                tail_cache[key] = t
            return tail_cache[key]
        def reject_prob(n1: int, n2: int, r1: int, r: int, p: float) -> float:
            pmf1 = get_pmf(n1, p)
            tail = get_tail(n2, p)
            prob = 0.0
            for x1 in range(r1 + 1, n1 + 1):
                needed = r + 1 - x1
                if needed <= 0:
                    # 第一阶段的有效数已经足够拒绝原假设，第二阶段概率为 1.0
                    prob += pmf1[x1] * 1.0
                elif needed > n2:
                    # 第二阶段需要的有效数超过了第二阶段的总人数，这不可能发生，概率为 0
                    prob += pmf1[x1] * 0.0
                else:
                    # 正常查表
                    prob += pmf1[x1] * tail[needed]
            return prob
        best_opt, best_min = None, None
        for n in range(2, n_max + 1):
            for n1 in range(1, n):
                if best_opt is not None and n1 >= best_opt["en0"]:
                    continue
                n2 = n - n1
                for r1 in range(0, n1):
                    pet0 = float(binom.cdf(r1, n1, p0))
                    en0 = n1 + (1.0 - pet0) * n2
                    if best_opt is not None and en0 >= best_opt["en0"]:
                        continue
                    for r in range(r1 + 1, n):
                        type1 = reject_prob(n1, n2, r1, r, p0)
                        if type1 > alpha:
                            continue
                        actual_pw = reject_prob(n1, n2, r1, r, p1)
                        if actual_pw < power:
                            continue
                        cand = {
                            "n1": n1, "n2": n2, "n": n,
                            "r1": r1, "r": r,
                            "pet0": pet0, "en0": en0,
                            "actual_alpha": type1, "actual_power": actual_pw,
                        }
                        if (best_opt is None
                            or cand["en0"] < best_opt["en0"]
                            or (abs(cand["en0"] - best_opt["en0"]) < 1e-12 and cand["n"] < best_opt["n"])):
                            best_opt = cand
                        if (best_min is None
                            or cand["n"] < best_min["n"]
                            or (cand["n"] == best_min["n"] and cand["en0"] < best_min["en0"])):
                            best_min = cand
                        break

        if best_opt is None or best_min is None:
            raise ValueError("在当前搜索范围内未找到满足约束的 Simon 二阶段设计。")

        selected = dict(best_opt if design_type == "optimal" else best_min)
        selected["design_type"] = design_type
        selected["total_sample_size_with_dropout"] = int(math.ceil(selected["n"] / (1 - dropout_rate)))
        selected["n1_with_dropout"] = int(math.ceil(selected["n1"] / (1 - dropout_rate)))
        selected["n2_with_dropout"] = int(math.ceil(selected["n2"] / (1 - dropout_rate)))
        selected["optimal_design"] = {
            "n1": best_opt["n1"], "n2": best_opt["n2"], "n": best_opt["n"],
            "r1": best_opt["r1"], "r": best_opt["r"],
            "en0": best_opt["en0"], "pet0": best_opt["pet0"],
            "actual_alpha": best_opt["actual_alpha"], "actual_power": best_opt["actual_power"],
        }
        selected["minimax_design"] = {
            "n1": best_min["n1"], "n2": best_min["n2"], "n": best_min["n"],
            "r1": best_min["r1"], "r": best_min["r"],
            "en0": best_min["en0"], "pet0": best_min["pet0"],
            "actual_alpha": best_min["actual_alpha"], "actual_power": best_min["actual_power"],
        }
        return selected

    # ---- 配对两组率 ----------------------------------------------------------
    def _paired_nam_n(
        self, p10: float, p01: float,
        alpha: float, power: float, test_type: str, margin: float,
    ) -> Tuple[float, float, float, float, float]:
        """Nam (1997) 型配对率样本量公式。

        Returns: (n_raw, p_d, psi, effect, inner)
        """
        p_d = p10 + p01
        psi = p10 - p01
        inner = max(p_d - psi ** 2, 0.0)

        z_alpha, z_beta = self._z_alpha_beta(alpha, power, test_type)

        if test_type == "差异性检验":
            if abs(psi) < 1e-12:
                raise ValueError("π₁₀ 与 π₀₁ 相同（率差 ψ=0），无法计算差异性样本量。")
            effect = abs(psi)
        elif test_type == "优效性检验":
            effect = psi - margin
            if effect <= 0:
                raise ValueError("优效性要求 π₁₀ − π₀₁ > Δ，请提高预期差异或缩小界值。")
        elif test_type == "非劣效性检验":
            effect = psi + margin
            if effect <= 0:
                raise ValueError("非劣效性要求 π₁₀ − π₀₁ > −Δ，请检查预期与界值。")
        elif test_type == "等效性检验":
            effect = margin - abs(psi)
            if effect <= 0:
                raise ValueError("等效性要求 |π₁₀ − π₀₁| < Δ，请调整预期或界值。")
        else:
            raise ValueError(f"不支持的检验类型: {test_type}")

        n_raw = (z_alpha * math.sqrt(p_d) + z_beta * math.sqrt(inner)) ** 2 / (effect ** 2)
        return n_raw, p_d, psi, effect, inner

    def _exact_mcnemar_diff_sample_size(
        self, p10: float, p01: float,
        alpha: float, power: float, dropout_rate: float,
    ) -> dict:
        """差异性：不一致对子集上条件二项精确检验（Exact McNemar）枚举搜索。"""
        p_d = p10 + p01
        p_alt = p10 / p_d
        if p_d <= 0 or p_alt <= 0 or p_alt >= 1:
            raise ValueError("Exact McNemar 需要有效的 π₁₀、π₀₁。")

        def power_at_n(n: int) -> float:
            d = max(1, int(round(n * p_d)))
            reject = [k for k in range(d + 1)
                      if float(binomtest(k, d, 0.5, alternative="two-sided").pvalue) <= alpha + 1e-15]
            return float(sum(binom.pmf(k, d, p_alt) for k in reject))

        n0 = int(max(5, math.ceil(self._paired_nam_n(p10, p01, alpha, power, "差异性检验", 0.0)[0])))
        for n in range(max(3, n0 - 20), 50001):
            if power_at_n(n) >= power - 1e-6:
                n_drop = int(math.ceil(n / (1 - dropout_rate)))
                return {
                    "n_pairs": n, "n_subjects": n,
                    "total_sample_size": n,
                    "n_pairs_with_dropout": n_drop,
                    "n_subjects_with_dropout": n_drop,
                    "total_sample_size_with_dropout": n_drop,
                    "method": "配对两组率：Exact McNemar 检验（条件二项精确）",
                    "method_ref": (
                        "Breslow NE, Day NE. Statistical Methods in Cancer Research, "
                        "Vol I (IARC, 1980), Section 4.4."
                    ),
                    "rate_note": (
                        "基于不一致对子集上的双侧二项精确检验，适用于 discordant pair 极少"
                        "或极小样本场景。大样本时与 McNemar 正态近似趋于一致。"
                    ),
                    "paired_psi": float(p10 - p01),
                    "paired_pd": float(p_d),
                }
        raise ValueError("Exact McNemar 在可搜索范围内未找到满足效能的样本量。")

    def _exact_conditional_margin_sample_size(
        self, p10: float, p01: float,
        alpha: float, power: float, dropout_rate: float,
        test_type: str, margin: float,
    ) -> dict:
        """优效/非劣：不一致对子集上单侧二项精确检验（敏感性/补充分析）。"""
        if test_type not in ("优效性检验", "非劣效性检验"):
            raise ValueError("exact conditional paired 仅适用于优效或非劣效。")
        if margin <= 0:
            raise ValueError("优效/非劣效须提供正界值 Δ。")
        p_d = p10 + p01
        p_alt = p10 / p_d
        if p_d <= 0 or p_alt <= 0 or p_alt >= 1:
            raise ValueError("无效的 π₁₀、π₀₁。")
        if test_type == "优效性检验":
            if p_d <= margin + 1e-12:
                raise ValueError("需 pD > Δ，以保证 H₀ 边界下 π₀₁ 非负。")
            p0_h0 = (p_d + margin) / (2.0 * p_d)
        else:
            if p_d <= margin + 1e-12:
                raise ValueError("需 pD > Δ，以保证 H₀ 边界有解。")
            p0_h0 = (p_d - margin) / (2.0 * p_d)
        if not (0 < p0_h0 < 1):
            raise ValueError("H₀ 边界下条件概率超出 (0,1)，请调整 pD 与 Δ。")

        def power_at_n(n: int) -> float:
            d = max(1, int(round(n * p_d)))
            reject = [k for k in range(d + 1)
                      if float(binomtest(k, d, p0_h0, alternative="greater").pvalue) <= alpha + 1e-15]
            return float(sum(binom.pmf(k, d, p_alt) for k in reject))

        n0 = int(max(5, math.ceil(
            self._paired_nam_n(p10, p01, alpha, power, test_type, margin)[0]
        )))
        for n in range(max(3, n0 - 40), 50001):
            if power_at_n(n) >= power - 1e-6:
                n_drop = int(math.ceil(n / (1 - dropout_rate)))
                return {
                    "n_pairs": n, "n_subjects": n,
                    "total_sample_size": n,
                    "n_pairs_with_dropout": n_drop,
                    "n_subjects_with_dropout": n_drop,
                    "total_sample_size_with_dropout": n_drop,
                    "n_treatment": n, "n_control": 0,
                    "n_treatment_with_dropout": n_drop,
                    "n_control_with_dropout": 0,
                    "method": "配对两组率：Exact Conditional Paired（条件二项单侧精确）",
                    "method_ref": (
                        "Nam J. Power and sample size for non-inferiority in paired "
                        "binary data. Stat Med, 1997."
                    ),
                    "rate_note": (
                        "适用于极小样本或敏感性/补充分析；"
                        "确证性主分析建议以 Nam Score 为主并交叉核对。"
                    ),
                    "paired_psi": float(p10 - p01),
                    "paired_pd": float(p_d),
                    "paired_margin": float(margin),
                }
        raise ValueError("exact conditional paired 在可搜索范围内未找到满足效能的样本量。")

    def compute_paired_proportion(
        self,
        p_discordant_treat: float, p_discordant_ctrl: float,
        alpha: float, power: float, dropout_rate: float,
        test_type: str,
        margin: float = 0.0,
        method: Optional[str] = None,
    ) -> dict:
        """配对两组率样本量。

        π₁₀、π₀₁：不一致配对比例；margin：优效/非劣/等效的率差界值。

        method:
          - None / "auto"：差异性自动选择 Nam/McNemar 或 Exact McNemar；其他默认 Nam Score
          - "nam_score"：Nam (1997) 正态近似（差异性即为 McNemar 正态近似）
          - "exact_mcnemar"：仅差异性（条件二项精确）
          - "exact_conditional_paired"：仅优效/非劣（条件二项精确）
        """
        p10, p01 = p_discordant_treat, p_discordant_ctrl
        if p10 <= 0 or p01 <= 0:
            raise ValueError("不一致配对比例（π₁₀、π₀₁）均须大于 0。")
        if p10 + p01 >= 1:
            raise ValueError("不一致配对比例之和须小于 1。")

        m = method or "auto"

        if m == "exact_conditional_paired":
            return self._exact_conditional_margin_sample_size(
                p10, p01, alpha, power, dropout_rate, test_type, margin,
            )
        if m == "exact_mcnemar" and test_type != "差异性检验":
            raise ValueError("Exact McNemar 仅适用于差异性检验（配对率 H₀: ψ=0）。")

        # 差异性自动选择
        if (m == "auto" or m == "exact_mcnemar") and test_type == "差异性检验":
            if m == "exact_mcnemar":
                return self._exact_mcnemar_diff_sample_size(
                    p10, p01, alpha, power, dropout_rate,
                )
            n_raw, p_d, *_ = self._paired_nam_n(p10, p01, alpha, power, "差异性检验", 0.0)
            n_est = int(math.ceil(n_raw))
            if n_est * p_d >= 10.0 and n_est >= 25:
                m_use = "nam_score"
            else:
                return self._exact_mcnemar_diff_sample_size(
                    p10, p01, alpha, power, dropout_rate,
                )
        else:
            m_use = "nam_score"

        # Nam Score 计算
        n_raw, p_d, psi, eff, inner = self._paired_nam_n(
            p10, p01, alpha, power, test_type, margin,
        )
        inner_note = f"ψ = π₁₀ − π₀₁ = {psi:.4f}，pD = π₁₀ + π₀₁ = {p_d:.4f}"

        if test_type == "等效性检验":
            method_label = "配对两组率：Nam Score TOST（双单侧正态近似，基于 discordant pairs）"
            method_ref = "Nam J. Stat Med, 1997."
            rate_note = (
                "等效性检验（TOST）为双单侧检验，每次单侧使用显著性水平 α。"
                "基于 Nam (1997) 提出的一致配对（discordant pair）方差结构，"
                "与独立样本的 Farrington–Manning TOST 思路一致。"
            )
        elif test_type == "差异性检验":
            method_label = "配对两组率：McNemar 检验（正态近似，基于 discordant pairs）"
            method_ref = "McNemar Q. Psychometrika, 1947; Nam J. Stat Med, 1997."
            rate_note = (
                f"配对二分类经典方法（基于不一致对子集）；{inner_note}。"
                "小样本或 discordant 很少（n·pD < 10）时请选 Exact McNemar。"
            )
        else:
            test_cn = {"优效性检验": "优效性", "非劣效性检验": "非劣效性"}.get(test_type, test_type)
            method_label = f"配对两组率：Nam Score 检验（{test_cn}，单侧，基于 discordant pairs）"
            method_ref = "Nam J. Stat Med, 1997."
            rate_note = (
                f"{test_cn}检验（单侧 α={alpha}），界值 Δ={margin:.4f}。{inner_note}。"
                "推荐以 Nam Score 为主方法；极小样本可用 exact conditional paired 做敏感性分析。"
            )

        n_pairs = max(1, int(math.ceil(n_raw)))
        n_drop = int(math.ceil(n_pairs / (1 - dropout_rate)))
        out = {
            "n_pairs": n_pairs,
            "n_subjects": n_pairs,
            "total_sample_size": n_pairs,
            "n_pairs_with_dropout": n_drop,
            "n_subjects_with_dropout": n_drop,
            "total_sample_size_with_dropout": n_drop,
            "n_treatment": n_pairs,
            "n_control": 0,
            "n_treatment_with_dropout": n_drop,
            "n_control_with_dropout": 0,
            "method": method_label,
            "method_ref": method_ref,
            "rate_note": rate_note,
            "paired_psi": float(psi),
            "paired_pd": float(p_d),
            "paired_margin": float(margin),
        }
        return out

    def compute_paired_mcnemar(
        self,
        p_discordant_treat: float, p_discordant_ctrl: float,
        alpha: float, power: float, dropout_rate: float,
        test_type: str,
    ) -> dict:
        """兼容旧接口：无界值 Nam 公式。"""
        return self.compute_paired_proportion(
            p_discordant_treat, p_discordant_ctrl,
            alpha, power, dropout_rate, test_type,
            margin=0.0, method="nam_score",
        )

    # ---- 三组率 --------------------------------------------------------------
    @staticmethod
    def _three_arm_ncp_pearson(ps: list) -> float:
        k = len(ps)
        p_bar = sum(ps) / k
        return sum((pi - p_bar) ** 2 for pi in ps) / (p_bar * (1 - p_bar))

    @staticmethod
    def _three_arm_min_expected_cell(ps: list, n: int) -> float:
        k = len(ps)
        p_bar = sum(ps) / k
        return min(n * p_bar, n * (1 - p_bar))

    def _three_arm_solve_n_ncx2(
        self, ncp_per_n: float, alpha: float, power: float, df: int
    ) -> int:
        if ncp_per_n < 1e-12:
            raise ValueError("三组有效率相同或效应量过小，无法计算样本量。")
        chi2_crit = float(chi2.ppf(1 - alpha, df))
        n = 5
        while n < 500000:
            if 1.0 - ncx2.cdf(chi2_crit, df, n * ncp_per_n) >= power:
                return n
            n += 1
        raise ValueError("在合理范围内未能找到满足效能的样本量，请检查参数。")






    def compute_three_arm_cochran_armitage(
        self,
        p1: float, p2: float, p3: float,
        alpha: float, power: float, dropout_rate: float,
        scores: tuple = (0.0, 1.0, 2.0),
    ) -> dict:
        """三组率 Cochran–Armitage 线性趋势检验（双侧，等额分配，非中心 χ² 法）。"""
        ps = [p1, p2, p3]
        k = 3
        sc = np.array(scores, dtype=float)
        if len(sc) != k:
            raise ValueError("趋势分数须为 3 个实数。")
        if not (sc[0] < sc[1] < sc[2]):
            raise ValueError("趋势分数须严格递增（如 0→1→2）。")
        p_bar = sum(ps) / k
        if p_bar <= 1e-9 or p_bar >= 1 - 1e-9:
            raise ValueError("三组平均有效率过于极端，无法计算。")
        t_bar = float(np.sum(sc) / k)
        trend_contrast = sum((sc[i] - t_bar) * (ps[i] - p_bar) for i in range(k))
        if abs(trend_contrast) < 1e-12:
            raise ValueError("在假定分数下，三组率相对均值的线性趋势分量为 0。")

        # CA 趋势检验的非中心参数（1 df）
        num = sum((sc[i] - t_bar) * (ps[i] - p_bar) for i in range(k)) ** 2
        denom = p_bar * (1.0 - p_bar) * sum((sc[i] - t_bar) ** 2 for i in range(k))
        ncp_per_n = num / denom if denom > 1e-15 else 0.0

        n = self._three_arm_solve_n_ncx2(ncp_per_n, alpha, power, df=1)
        n_per_group_drop = int(math.ceil(n / (1 - dropout_rate)))
        min_cell = self._three_arm_min_expected_cell(ps, n)
        return {
            "n_per_group": n,
            "total_sample_size": k * n,
            "n_per_group_with_dropout": n_per_group_drop,
            "total_sample_size_with_dropout": k * n_per_group_drop,
            "method": "三组独立样本率：Cochran–Armitage 趋势检验（非中心 χ² 法，双侧）",
            "method_ref": "Cochran WG. Biometrics, 1954; Armitage P. Biometrics, 1955.",
            "rate_note": (
                f"假定三组有序且关注线性趋势；分数为 ({scores[0]:g}, {scores[1]:g}, {scores[2]:g})。"
                "基于渐近 Z 统计量的非中心 χ² 法估计功效。"
            ),
            "three_arm_method": "cochran_armitage",
            "min_expected_cell_approx": float(min_cell),
            "ca_scores": scores,
        }

    def compute_three_arm_proportion(
        self,
        p1: float, p2: float, p3: float,
        alpha: float, power: float, dropout_rate: float,
        method: Optional[str] = None,
    ) -> dict:
        """三组率差异性检验样本量（k×2 齐性，Pearson χ² 非中心法，等额分配）。"""
        ps = [p1, p2, p3]
        k, df = 3, 2
        p_bar = sum(ps) / k
        if p_bar <= 1e-9 or p_bar >= 1 - 1e-9:
            raise ValueError("三组平均有效率过于极端，无法计算。")
        ncp_p = self._three_arm_ncp_pearson(ps)
        if ncp_p < 1e-12:
            raise ValueError("三组有效率相同或效应量过小。")

        n = self._three_arm_solve_n_ncx2(ncp_p, alpha, power, df)

        n_per_group_drop = int(math.ceil(n / (1 - dropout_rate)))
        min_cell = self._three_arm_min_expected_cell(ps, n)
        normal_zone = all(0.20 <= pi <= 0.80 for pi in ps)
        if not normal_zone or min_cell < 5.0:
            rate_note = (
                "基于渐近 Pearson χ² 的非中心参数法。"
                "当前至少一组率不在 [0.2, 0.8] 或期望频数 < 5，"
                "建议在分析阶段使用确切检验（需 StatXact 等专业软件）进行复核。"
            )
        else:
            rate_note = "基于渐近 Pearson χ² 的非中心参数法；正态近似可靠。"
        return {
            "n_per_group": n,
            "total_sample_size": k * n,
            "n_per_group_with_dropout": n_per_group_drop,
            "total_sample_size_with_dropout": k * n_per_group_drop,
            "method": "三组独立样本率：Pearson χ² 齐性检验（非中心 χ² 法，等额分配）",
            "method_ref": "Chow SC, et al. (3rd Ed), Section 4.5.",
            "rate_note": rate_note,
            "three_arm_method": "pearson",
            "min_expected_cell_approx": float(min_cell),
            "ncp_per_n": float(ncp_p),
        }

    # ---- 通用格式化 ----------------------------------------------------------
    def _format_result(self, nt: float, nc: float, dropout: float) -> dict:
        nt_ceil = int(math.ceil(nt))
        nc_ceil = int(math.ceil(nc))
        nt_drop = int(math.ceil(nt_ceil / (1 - dropout)))
        nc_drop = int(math.ceil(nc_ceil / (1 - dropout)))
        return {
            "n_treatment": nt_ceil,
            "n_control": nc_ceil,
            "total_sample_size": nt_ceil + nc_ceil,
            "n_treatment_with_dropout": nt_drop,
            "n_control_with_dropout": nc_drop,
            "total_sample_size_with_dropout": nt_drop + nc_drop,
        }
