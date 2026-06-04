# calc_proportion.py
from dataclasses import dataclass
from typing import Any, Dict, Literal, Optional
import math
import numpy as np
from scipy.stats import norm, binom, chi2, ncx2, beta, fisher_exact, binomtest, chi2_contingency

TestType = Literal["差异性检验", "优效性检验", "非劣效性检验", "等效性检验"]
Direction = Literal["higher", "lower"]  # higher: 高值/高有效率更优；lower: 低事件率更优

@dataclass(frozen=True)
class ProportionInput:
    p_treatment: float
    p_control: Optional[float] = None  
    alpha: float = 0.05
    power: float = 0.8
    test_type: TestType = "优效性检验"
    allocation_ratio: float = 1.0
    dropout_rate: float = 0.2
    margin: float = 0.0
    direction: Direction = "higher"
    # 双臂独立两组：差异性 / 优非劣等 的算法选择（None=自动按率与 np 规则）
    two_arm_diff_method: Optional[str] = None
    two_arm_margin_method: Optional[str] = None

class ProportionCalculator:
    @staticmethod
    def _rates_in_normal_zone(pt: float, pc: float) -> bool:
        """正态近似 Z 检验常用适用区间：两组率均在 [0.2, 0.8]。"""
        return 0.2 <= pt <= 0.8 and 0.2 <= pc <= 0.8

    @staticmethod
    def _expected_counts_ok(pt: float, pc: float, r: float, n_t: float) -> bool:
        """基于合并方差公式得到的试验组样本量，检验四格期望频数是否均 ≥5（Cochran 常用规则）。"""
        if n_t <= 0:
            return False
        n_t = int(math.ceil(n_t))
        n_c = int(math.ceil(n_t * r))
        e = (n_t * pt, n_t * (1 - pt), n_c * pc, n_c * (1 - pc))
        return min(e) >= 5.0 - 1e-9

    @staticmethod
    def _clip_p(p: float) -> float:
        return min(max(p, 1e-6), 1.0 - 1e-6)

    def _v0_projection_superiority(self, pt: float, pc: float, margin: float, r: float) -> float:
        """H0 边界 p_T − p_C = margin 上，(p_T,p_C) 到 (pt,pc) 的投影方差（近似 FM 零假设方差）。"""
        x = (pt + pc - margin) / 2.0
        p_tilde_c = self._clip_p(x)
        p_tilde_t = self._clip_p(x + margin)
        return p_tilde_t * (1 - p_tilde_t) + (p_tilde_c * (1 - p_tilde_c)) / r

    def _v0_projection_noninferiority(self, pt: float, pc: float, margin: float, r: float) -> float:
        """H0 边界 p_T − p_C = −margin 上的投影方差。"""
        x = (pt + pc + margin) / 2.0
        p_tilde_c = self._clip_p(x)
        p_tilde_t = self._clip_p(x - margin)
        return p_tilde_t * (1 - p_tilde_t) + (p_tilde_c * (1 - p_tilde_c)) / r

    @staticmethod
    def _fisher_two_sample_power(n1: int, n2: int, p1: float, p2: float, alpha: float) -> float:
        """双侧 Fisher 精确检验（2×2）在 H0: p1=p2 下的渐近功效（枚举所有表）。"""
        s = 0.0
        for a in range(n1 + 1):
            pmf_a = float(binom.pmf(a, n1, p1))
            if pmf_a < 1e-15:
                continue
            b = n1 - a
            for c in range(n2 + 1):
                pmf_c = float(binom.pmf(c, n2, p2))
                if pmf_c < 1e-20:
                    continue
                d = n2 - c
                try:
                    pv = float(fisher_exact([[a, b], [c, d]], alternative="two-sided").pvalue)
                except ValueError:
                    continue
                if pv <= alpha:
                    s += pmf_a * pmf_c
        return float(s)

    def _fisher_exact_sample_size_diff(
        self,
        pt: float,
        pc: float,
        r: float,
        alpha: float,
        power: float,
        dropout: float,
    ) -> Dict[str, Any]:
        """差异性：Fisher 精确检验，枚举 n 搜索最小样本量（适用于中小样本）。"""
        z_a = float(norm.ppf(1 - alpha / 2.0))
        z_b = float(norm.ppf(power))
        effect = abs(pt - pc)
        if effect < 1e-12:
            raise ValueError("两组预期率相同，无需计算样本量。")
        p_avg = (pt + r * pc) / (1 + r)
        v0 = p_avg * (1 - p_avg) * (1 + 1 / r)
        v1 = pt * (1 - pt) + (pc * (1 - pc)) / r
        n0 = int(max(8, math.ceil(((z_a * math.sqrt(v0) + z_b * math.sqrt(v1)) ** 2) / (effect**2))))
        n_start = max(5, n0 - 30)
        n_max = 8000
        for n1 in range(n_start, n_max + 1):
            n2 = int(math.ceil(n1 * r))
            pw = self._fisher_two_sample_power(n1, n2, pt, pc, alpha)
            if pw >= power:
                out = self._format_result(float(n1), float(n2), dropout)
                out["method"] = "两独立样本率差：Fisher 精确检验（双侧，枚举功效）"
                out["rate_zone"] = "fisher"
                out["rate_note"] = (
                    "基于 Fisher 精确检验枚举计算样本量；大样本时与 Pearson 卡方/合并方差 Z 接近。"
                    "若 n 搜索触及上限，请改用合并方差 Z 或反正弦法。"
                )
                out["actual_power_at_alt"] = float(pw)
                return out
        raise ValueError(
            "Fisher 精确检验在可搜索范围内未找到满足效能的样本量；请放宽效应或改用正态近似法。"
        )

    def _v0_wald_simple_margin(
        self,
        pt: float,
        pc: float,
        margin: float,
        r: float,
        test: str,
        direction: str,
    ) -> float:
        """Wald 型：H₀ 边界上直接代入方差（非投影），用于大样本简化。"""
        if test == "优效性检验":
            if direction == "higher":
                delta = margin
            else:
                delta = -margin
        elif test == "非劣效性检验":
            delta = -margin
        elif test == "等效性检验":
            pt_sup = min(max(pc + margin, 0.001), 0.999)
            pt_inf = min(max(pc - margin, 0.001), 0.999)
            v_sup = pt_sup * (1 - pt_sup) + (pc * (1 - pc)) / r
            v_inf = pt_inf * (1 - pt_inf) + (pc * (1 - pc)) / r
            return float(max(v_sup, v_inf))
        else:
            raise ValueError("仅支持优效/非劣/等效的 Wald 边界方差。")

        pt0 = min(max(pc + delta, 0.001), 0.999)
        return float(pt0 * (1 - pt0) + (pc * (1 - pc)) / r)

    def compute_two_arm(self, params: ProportionInput) -> Dict[str, Any]:
        pt, pc, r, margin, test = (
            params.p_treatment,
            params.p_control,
            params.allocation_ratio,
            params.margin,
            params.test_type,
        )
        direction = params.direction
        if test == "差异性检验":
            delta, z_alpha, effect = 0.0, norm.ppf(1 - params.alpha / 2), abs(pt - pc)
        elif test == "优效性检验":
            # 高值更优: H1: pT - pC > Δ
            # 低值更优: H1: pC - pT > Δ  (例如不良事件发生率降低)
            if direction == "higher":
                delta, z_alpha, effect = margin, norm.ppf(1 - params.alpha), (pt - pc) - margin
            elif direction == "lower":
                delta, z_alpha, effect = -margin, norm.ppf(1 - params.alpha), (pc - pt) - margin
            else:
                raise ValueError(f"不支持的 direction: {direction}")
        elif test == "非劣效性检验":
            delta, z_alpha, effect = -margin, norm.ppf(1 - params.alpha), (pt - pc) + margin
        elif test == "等效性检验":
            delta, z_alpha, effect = margin, norm.ppf(1 - params.alpha), margin - abs(pt - pc)
        else:
            raise ValueError(f"不支持的检验类型: {test}")

        z_beta = norm.ppf(1 - (1 - params.power) / 2.0) if test == "等效性检验" else norm.ppf(params.power)
        if effect <= 0:
            raise ValueError(f"当前设定下效应量 ≤ 0，无法满足【{test}】要求。")

        normal_zone = self._rates_in_normal_zone(pt, pc)
        v1 = pt * (1 - pt) + (pc * (1 - pc)) / r

        # ---------- 差异性：合并方差 Z / Fisher 精确 / 反正弦 ----------
        if test == "差异性检验":
            dm = params.two_arm_diff_method or "auto"
            if dm == "auto":
                p_avg = (pt + r * pc) / (1 + r)
                v0 = p_avg * (1 - p_avg) * (1 + 1 / r)
                nt_plan = ((z_alpha * math.sqrt(v0) + z_beta * math.sqrt(v1)) ** 2) / (effect ** 2)
                # 按合并方差估计的样本量下四格期望频数均 ≥5 → Pooled Z（含罕见病但 n 足够大）；否则 Fisher
                if self._expected_counts_ok(pt, pc, r, nt_plan):
                    dm = "pooled_z"
                else:
                    dm = "fisher_exact"

            if dm == "fisher_exact":
                return self._fisher_exact_sample_size_diff(
                    pt, pc, r, params.alpha, params.power, params.dropout_rate
                )

            if dm == "arcsine":
                st = math.asin(math.sqrt(self._clip_p(pt)))
                sc = math.asin(math.sqrt(self._clip_p(pc)))
                d_arc = st - sc
                if abs(d_arc) < 1e-12:
                    raise ValueError("反正弦变换后效应量过小，无法计算样本量。")
                denom = d_arc ** 2
                nt = ((z_alpha * math.sqrt((1 + 1 / r) / 4.0) + z_beta * math.sqrt((1 + 1 / r) / 4.0)) ** 2) / denom
                out = self._format_result(nt, nt * r, params.dropout_rate)
                out["method"] = "两独立样本率差：反正弦变换（方差稳定化，极端率备选）"
                out["rate_zone"] = "extreme"
                out["rate_note"] = (
                    "任选反正弦法；适用于极端率。分析阶段可与 Pearson 卡方、Fisher 精确互为补充。"
                )
                return out

            # pooled_z
            p_avg = (pt + r * pc) / (1 + r)
            v0 = p_avg * (1 - p_avg) * (1 + 1 / r)
            method_label = "两独立样本率差：正态近似 Z 检验（合并方差，与 Pearson 卡方等价）"
            nt_raw = ((z_alpha * math.sqrt(v0) + z_beta * math.sqrt(v1)) ** 2) / (effect ** 2)
            nt_int = int(math.ceil(nt_raw))
            nc_int = int(math.ceil(nt_int * r))
            exp_ok = self._expected_counts_ok(pt, pc, r, nt_raw)
            if exp_ok and normal_zone:
                rate_note = (
                    "合并方差 Z（双侧）与 Pearson χ² 独立性检验等价；"
                    "按当前假设估计样本量下四格期望频数均 ≥5（Cochran 常用规则）。"
                )
                zone = "normal"
            elif exp_ok:
                rate_note = (
                    "合并方差 Z；至少一组率不在 [0.2,0.8]，但按当前假设估计样本量下期望频数仍 ≥5。"
                    "若确证性分析仍担心近似，可复核 Fisher 或反正弦。"
                )
                zone = "extreme_pooled"
            else:
                rate_note = (
                    "合并方差 Z（用户指定）；按当前假设估计样本量下部分期望频数 <5，"
                    "正态近似与 χ² 可能偏倚，分析阶段建议 Fisher 精确或改用上方「自动」策略。"
                )
                zone = "extreme_pooled_sparse"
            out = self._format_result(nt_raw, nt_raw * r, params.dropout_rate)
            out["method"] = method_label
            out["rate_zone"] = zone
            out["rate_note"] = rate_note
            out["expected_min_cell"] = float(
                min(nt_int * pt, nt_int * (1 - pt), nc_int * pc, nc_int * (1 - pc))
            )
            return out

        # ---------- 优效 / 非劣 / 等效 ----------
        mm = params.two_arm_margin_method or "auto"
        chan_note = ""
        chan_suffix = False
        if mm == "auto":
            mm = "fm_score"
            if not normal_zone:
                chan_note = (
                    "说明：极端率或极小样本下，完整 Chan 无条件精确检验需专用软件（如 StatXact/Cytel）；"
                    "此处样本量按 FM Score 投影方差近似，注册试验建议独立复核。"
                )
        elif mm == "unconditional_chan":
            mm = "fm_score"
            chan_note = (
                "说明：完整 Chan 无条件精确检验需专用软件（如 StatXact/Cytel）；"
                "此处样本量按 FM Score 投影方差近似，注册试验建议独立复核。"
            )
            chan_suffix = True

        if mm == "wald_unpooled":
            v0 = self._v0_wald_simple_margin(pt, pc, margin, r, test, direction)
            method_label = "两独立样本率：Wald 型非合并方差（H₀ 边界直接代入，大样本近似）"
            rate_note = (
                "Wald 型边界方差；大样本下与 FM Score 接近。率极偏或小样本优先 FM Score 或无条件精确类方法。"
            )
            nt = ((z_alpha * math.sqrt(v0) + z_beta * math.sqrt(v1)) ** 2) / (effect ** 2)
            out = self._format_result(nt, nt * r, params.dropout_rate)
            out["method"] = method_label
            out["rate_zone"] = "wald"
            out["rate_note"] = rate_note
            return out

        # fm_score（含 unconditional_chan 的数值路径）
        if normal_zone:
            pt0 = min(max(pc + delta, 0.001), 0.999)
            v0 = pt0 * (1 - pt0) + (pc * (1 - pc)) / r
            if test == "等效性检验":
                method_label = "两独立样本率：FM TOST（两个单侧 FM Score，H₀ 边界方差；ICH 常用正态近似框架）"
                rate_note = (
                    "等效性采用双单侧检验（TOST）；两组率均在 [0.2,0.8]，与 SAS/East 等 FM 思路一致。"
                )
            else:
                method_label = "两独立样本率：FM Score（H₀ 边界方差，近似 Farrington–Manning）"
                rate_note = "两组率均在 [0.2, 0.8]，药政申报常用基准之一（与 SAS/East 等实现思路一致）。"
        else:
            if test == "优效性检验":
                if direction == "higher":
                    v0 = self._v0_projection_superiority(pt, pc, margin, r)
                else:
                    v0 = self._v0_projection_noninferiority(pt, pc, margin, r)
                method_label = "两独立样本率：FM Score（极端率投影方差）"
            elif test == "非劣效性检验":
                v0 = self._v0_projection_noninferiority(pt, pc, margin, r)
                method_label = "两独立样本率：FM Score（极端率投影方差）"
            else:
                v0 = max(
                    self._v0_projection_superiority(pt, pc, margin, r),
                    self._v0_projection_noninferiority(pt, pc, margin, r),
                )
                method_label = "两独立样本率：FM TOST（极端率：两侧 H₀ 边界投影方差取较大者，保守）"
            rate_note = (
                "至少一组率不在 [0.2, 0.8]，已采用 FM Score 投影方差；小样本或确证性试验可复核无条件精确方法。"
            )
        if chan_note:
            rate_note = chan_note + " " + rate_note
        if chan_suffix:
            method_label = method_label + "（Chan 场景参考）"

        nt = ((z_alpha * math.sqrt(v0) + z_beta * math.sqrt(v1)) ** 2) / (effect ** 2)
        out = self._format_result(nt, nt * r, params.dropout_rate)
        out["method"] = method_label
        out["rate_zone"] = "normal" if normal_zone else "extreme"
        out["rate_note"] = rate_note
        return out

    @staticmethod
    def _cp_ci(x: int, n: int, alpha: float) -> tuple[float, float]:
        """Clopper-Pearson exact CI for a single proportion."""
        if n <= 0:
            raise ValueError("n 必须为正整数。")
        if not (0 <= x <= n):
            raise ValueError("x 必须位于 [0,n]。")
        if x == 0:
            lo = 0.0
        else:
            lo = float(beta.ppf(alpha / 2.0, x, n - x + 1))
        if x == n:
            hi = 1.0
        else:
            hi = float(beta.ppf(1 - alpha / 2.0, x + 1, n - x))
        return lo, hi

    @staticmethod
    def _wilson_ci(x: int, n: int, alpha: float) -> tuple[float, float]:
        """Wilson score CI for a single proportion."""
        if n <= 0:
            raise ValueError("n 必须为正整数。")
        if not (0 <= x <= n):
            raise ValueError("x 必须位于 [0,n]。")
        z = float(norm.ppf(1 - alpha / 2.0))
        phat = x / n
        denom = 1.0 + (z**2) / n
        center = (phat + (z**2) / (2 * n)) / denom
        half = (z / denom) * math.sqrt((phat * (1 - phat) + (z**2) / (4 * n)) / n)
        lo = max(0.0, center - half)
        hi = min(1.0, center + half)
        return lo, hi

    def compute_single_arm(self, params: ProportionInput, method: Literal["score", "exact"] = "score") -> Dict[str, Any]:
        """单臂率样本量（功效分析）"""
        pt, p0, margin, test = params.p_treatment, params.p_control, params.margin, params.test_type
        if p0 is None:
            raise ValueError("单臂设计需要提供历史/目标有效率 P0。")
        if not (0 < pt < 1 and 0 < p0 < 1):
            raise ValueError("P0/P1 必须位于 (0,1)。")
        if not (0 < params.alpha < 1 and 0 < params.power < 1):
            raise ValueError("alpha/power 必须位于 (0,1)。")

        # 统一定义 H0 边界点 p0_bound（用于 z_alpha 的方差）
        if test == "差异性检验":
            alpha_side = "two-sided"
            p0_bound = p0
            z_alpha = norm.ppf(1 - params.alpha / 2.0)
            effect = abs(pt - p0)
        elif test == "优效性检验":
            alpha_side = "one-sided"
            p0_bound = p0 + margin
            z_alpha = norm.ppf(1 - params.alpha)
            effect = pt - (p0 + margin)
        elif test == "非劣效性检验":
            alpha_side = "one-sided"
            p0_bound = p0 - margin
            z_alpha = norm.ppf(1 - params.alpha)
            effect = pt - (p0 - margin)
        elif test == "等效性检验":
            # TOST：alpha 为单侧水平（常用 0.05），此处保持与现有实现一致
            alpha_side = "tost"
            p0_bound = p0  # 精度上等效两侧边界不同；近似下仍用 p0 做 zα 方差
            z_alpha = norm.ppf(1 - params.alpha)
            effect = margin - abs(pt - p0)
        else:
            raise ValueError(f"不支持的检验类型: {test}")

        if test == "等效性检验":
            z_beta = norm.ppf(1 - (1 - params.power) / 2.0)
        else:
            z_beta = norm.ppf(params.power)

        if effect <= 0:
            raise ValueError(f"当前设定下效应量 ≤ 0，无法满足【{test}】要求。")

        p0_bound = min(max(p0_bound, 0.001), 0.999)

        if method == "score":
            # Score 正态近似：zα 用 H0 方差，zβ 用 H1 方差
            n = ((z_alpha * math.sqrt(p0_bound * (1 - p0_bound)) + z_beta * math.sqrt(pt * (1 - pt))) ** 2) / (effect ** 2)
            out = self._format_result(n, 0, params.dropout_rate)
            out["method"] = "单臂率：Score 正态近似（zα 使用 H0 方差，zβ 使用 H1 方差）"
            out["alpha_side"] = alpha_side
            out["p0_bound"] = float(p0_bound)
            return out

        if method != "exact":
            raise ValueError("method 仅支持 'score' 或 'exact'。")

        # 精确二项：通过枚举 x=0..n 的精确检验规则计算 power，并以近似 n0 作为起点减少搜索区间。
        n0 = int(self._format_result(((z_alpha * math.sqrt(p0_bound * (1 - p0_bound)) + z_beta * math.sqrt(pt * (1 - pt))) ** 2) / (effect ** 2), 0, 0)["total_sample_size"])
        n_start = max(5, n0 - 50)
        n_max = 5000

        # 缓存 binom.sf / cdf，避免重复计算
        sf_cache: dict[tuple[int, float, int], float] = {}
        cdf_cache: dict[tuple[int, float, int], float] = {}

        def sf(n: int, p: float, k: int) -> float:
            key = (n, p, k)
            if key not in sf_cache:
                sf_cache[key] = float(binom.sf(k - 1, n, p))
            return sf_cache[key]

        def cdf(n: int, p: float, k: int) -> float:
            key = (n, p, k)
            if key not in cdf_cache:
                cdf_cache[key] = float(binom.cdf(k, n, p))
            return cdf_cache[key]

        def exact_p_two_sided(n: int, p: float, x: int) -> float:
            # 常用的“倍增尾概率”双侧精确二项 p 值（与前面讨论保持一致）
            phat = x / n
            if phat >= p:
                pv = 2.0 * sf(n, p, x)
            else:
                pv = 2.0 * cdf(n, p, x)
            return min(pv, 1.0)

        def exact_p_one_sided_greater(n: int, p: float, x: int) -> float:
            return sf(n, p, x)

        def exact_p_one_sided_less(n: int, p: float, x: int) -> float:
            return cdf(n, p, x)

        def power_from_rejection_mask(n: int, p_h1: float, reject: list[bool]) -> float:
            # power = Σ I(reject[x]) * BinomPMF(x|n,p1)
            pmf = [float(binom.pmf(x, n, p_h1)) for x in range(n + 1)]
            return float(sum(pmf[x] for x in range(n + 1) if reject[x]))

        def alpha_from_rejection_mask(n: int, p_h0: float, reject: list[bool]) -> float:
            pmf = [float(binom.pmf(x, n, p_h0)) for x in range(n + 1)]
            return float(sum(pmf[x] for x in range(n + 1) if reject[x]))

        for n in range(n_start, n_max + 1):
            reject = [False] * (n + 1)

            if test == "差异性检验":
                for x in range(n + 1):
                    if exact_p_two_sided(n, p0_bound, x) <= params.alpha:
                        reject[x] = True
            elif test == "优效性检验":
                # H0: p <= p0_bound，右尾拒绝
                for x in range(n + 1):
                    if exact_p_one_sided_greater(n, p0_bound, x) <= params.alpha:
                        reject[x] = True
            elif test == "非劣效性检验":
                # H0: p <= p0_bound，右尾拒绝（非劣：p > p0-Δ）
                for x in range(n + 1):
                    if exact_p_one_sided_greater(n, p0_bound, x) <= params.alpha:
                        reject[x] = True
            elif test == "等效性检验":
                # TOST：同时通过两侧单侧检验
                p_lo = min(max(p0 - margin, 0.001), 0.999)
                p_hi = min(max(p0 + margin, 0.001), 0.999)
                for x in range(n + 1):
                    p1v = exact_p_one_sided_greater(n, p_lo, x)  # 检验 p > p_lo
                    p2v = exact_p_one_sided_less(n, p_hi, x)     # 检验 p < p_hi
                    if p1v <= params.alpha and p2v <= params.alpha:
                        reject[x] = True
            else:
                raise ValueError(f"不支持的检验类型: {test}")

            pw = power_from_rejection_mask(n, pt, reject)
            if pw >= params.power:
                out = self._format_result(n, 0, params.dropout_rate)
                out["method"] = "单臂率：精确二项（Exact）"
                out["alpha_side"] = alpha_side
                out["p0_bound"] = float(p0_bound)
                out["actual_power_at_p1"] = float(pw)
                out["actual_alpha_at_p0"] = float(alpha_from_rejection_mask(n, p0_bound, reject))
                out["n0_start"] = int(n0)
                out["n_start"] = int(n_start)
                return out

        raise ValueError("精确二项在搜索上限内未找到满足效能的样本量，请检查参数或改用近似法。")

    def compute_single_arm_precision(
        self,
        p: float,
        ci_level: float,
        half_width: float,
        method: Literal["clopper_pearson", "wilson", "wald"] = "clopper_pearson",
        strategy: Literal["point"] = "point",
        n_max: int = 200000,
    ) -> Dict[str, Any]:
        """单臂率精度分析：给定目标 CI 半宽，反推最小样本量。"""
        if not (0 < p < 1):
            raise ValueError("p 必须位于 (0,1)。")
        if not (0 < ci_level < 1):
            raise ValueError("CI 置信水平必须位于 (0,1)。")
        if not (0 < half_width < 0.5):
            raise ValueError("half_width 必须位于 (0,0.5)。")
        alpha = 1.0 - ci_level

        if method == "wald":
            z = float(norm.ppf(1 - alpha / 2.0))
            n = max(1, int(math.ceil((z**2) * p * (1 - p) / (half_width**2))))
            hw = z * math.sqrt(p * (1 - p) / n)
            x = int(round(n * p))
            return {
                "n": n,
                "x_reference": x,
                "p_reference": float(p),
                "ci_level": float(ci_level),
                "half_width_target": float(half_width),
                "half_width_achieved": float(hw),
                "ci_lower": float(p - hw),
                "ci_upper": float(p + hw),
                "method": "单臂率精度分析：正态近似法 (Wald)",
            }

        # 用 Wald 半宽粗估作为起点（避免从 1 枚举）
        z = float(norm.ppf(1 - alpha / 2.0))
        n0 = int(math.ceil((z**2) * p * (1 - p) / (half_width**2)))
        n = max(5, n0 - 50)

        def interval(n_: int) -> tuple[float, float, int]:
            x = int(round(n_ * p)) if strategy == "point" else int(round(n_ * p))
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
                return {
                    "n": n,
                    "x_reference": x,
                    "p_reference": float(p),
                    "ci_level": float(ci_level),
                    "half_width_target": float(half_width),
                    "half_width_achieved": float(hw),
                    "ci_lower": float(lo),
                    "ci_upper": float(hi),
                    "method": "单臂率精度分析：Clopper-Pearson 精确CI" if method == "clopper_pearson" else "单臂率精度分析：Wilson Score CI",
                }
            n += 1

        raise ValueError("在搜索上限内未找到满足精度要求的样本量，请检查参数。")

    def compute_two_arm_precision(
        self,
        p_t: float,
        p_c: float,
        r: float,
        ci_level: float,
        half_width: float,
        method: Literal["wald", "newcombe_wilson"] = "wald",
        n_max: int = 250000,
    ) -> Dict[str, Any]:
        """两臂率差 (p_T−p_C) 的置信区间半宽不超过给定值时的最小样本量（试验组 n_T，对照组 n_C=r·n_T）。"""
        if not (0 < p_t < 1 and 0 < p_c < 1):
            raise ValueError("两组预期率须位于 (0,1)。")
        if not (0 < ci_level < 1):
            raise ValueError("置信水平须位于 (0,1)。")
        if not (0 < half_width < 1):
            raise ValueError("率差半宽须介于 (0,1)。")
        alpha = 1.0 - ci_level
        z = float(norm.ppf(1 - alpha / 2.0))

        if method == "wald":
            n_t = max(
                1,
                int(
                    math.ceil(
                        (z**2) * (p_t * (1 - p_t) + p_c * (1 - p_c) / r) / (half_width**2)
                    )
                ),
            )
            n_c = int(math.ceil(n_t * r))
            se = math.sqrt(p_t * (1 - p_t) / n_t + p_c * (1 - p_c) / n_c)
            hw = z * se
            return {
                "n_treatment": n_t,
                "n_control": n_c,
                "total_sample_size": n_t + n_c,
                "ci_level": float(ci_level),
                "half_width_target": float(half_width),
                "half_width_achieved": float(hw),
                "method": "两臂率差精度：Wald 近似（正态近似率差 SE）",
                "p_treatment": float(p_t),
                "p_control": float(p_c),
            }

        # Newcombe (1998) method 10：基于两侧 Wilson 区间构造率差区间，取保守最大单侧偏离为精度准则
        for n_t in range(5, n_max + 1):
            n_c = int(math.ceil(n_t * r))
            x_t = int(round(n_t * p_t))
            x_c = int(round(n_c * p_c))
            x_t = min(max(x_t, 0), n_t)
            x_c = min(max(x_c, 0), n_c)
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
                    "n_treatment": n_t,
                    "n_control": n_c,
                    "total_sample_size": n_t + n_c,
                    "ci_level": float(ci_level),
                    "half_width_target": float(half_width),
                    "half_width_achieved": float(asym),
                    "ci_lower_diff": float(L),
                    "ci_upper_diff": float(U),
                    "method": "两臂率差精度：Newcombe-Wilson Score（方法 10）",
                    "p_treatment": float(p_t),
                    "p_control": float(p_c),
                }

        raise ValueError("在搜索上限内未找到满足精度要求的样本量，请检查参数或改用 Wald。")

    def compute_simon_two_stage(
        self,
        p0: float,
        p1: float,
        alpha: float,
        power: float,
        dropout_rate: float,
        design_type: Literal["optimal", "minimax"] = "optimal",
        n_max: int = 100,
    ) -> dict:
        """
        精确 Simon 两阶段单臂设计搜索。

        约束:
          - Type I error <= alpha  (在 p0 下)
          - Power >= power         (在 p1 下)

        design_type:
          - "optimal":   在可行设计中最小化 EN(p0)
          - "minimax":   在可行设计中最小化总样本量 n
        """
        if not (0 < p0 < 1 and 0 < p1 < 1):
            raise ValueError("P0/P1 必须位于 (0,1)。")
        if p1 <= p0:
            raise ValueError("目标有效率 P1 必须大于无效界值 P0。")
        if not (0 < alpha < 1 and 0 < power < 1):
            raise ValueError("alpha/power 必须位于 (0,1)。")
        if not (0 <= dropout_rate < 1):
            raise ValueError("脱落率必须位于 [0,1)。")

        if design_type not in ("optimal", "minimax"):
            raise ValueError("design_type 仅支持 'optimal' 或 'minimax'。")

        pmf_cache = {}
        tail_cache = {}

        def get_pmf(n1: int, p: float):
            key = (n1, p)
            if key not in pmf_cache:
                pmf_cache[key] = [float(binom.pmf(x, n1, p)) for x in range(n1 + 1)]
            return pmf_cache[key]

        def get_tail(n2: int, p: float):
            # tail[k] = P(X >= k), k in [0, n2+1]
            key = (n2, p)
            if key not in tail_cache:
                tail = [1.0] * (n2 + 2)
                for k in range(1, n2 + 1):
                    tail[k] = float(binom.sf(k - 1, n2, p))
                tail[n2 + 1] = 0.0
                tail_cache[key] = tail
            return tail_cache[key]

        def reject_prob(n1: int, n2: int, r1: int, r: int, p: float) -> float:
            """P(最终拒绝H0 | p)"""
            pmf1 = get_pmf(n1, p)
            tail = get_tail(n2, p)
            prob = 0.0
            for x1 in range(r1 + 1, n1 + 1):
                # 最终拒绝条件: x1 + x2 > r  -> x2 >= r + 1 - x1
                needed = r + 1 - x1
                if needed <= 0:
                    tail_prob = 1.0
                elif needed > n2:
                    tail_prob = 0.0
                else:
                    tail_prob = tail[needed]
                prob += pmf1[x1] * tail_prob
            return prob

        best_opt = None
        best_min = None

        # n 为两阶段总理论样本量（不含脱落）
        for n in range(2, n_max + 1):
            for n1 in range(1, n):
                # 对 optimal 的安全剪枝：EN0 >= n1
                if best_opt is not None and n1 >= best_opt["en0"]:
                    continue

                n2 = n - n1
                for r1 in range(0, n1):
                    pet0 = float(binom.cdf(r1, n1, p0))  # 第一阶段在 p0 下提前终止概率
                    en0 = n1 + (1.0 - pet0) * n2

                    if best_opt is not None and en0 >= best_opt["en0"]:
                        continue

                    # r 为最终拒绝阈值：总有效例数 > r 时拒绝 H0
                    # r 必须 >= r1 + 1，否则第一阶段继续后几乎必拒绝不合理
                    for r in range(r1 + 1, n):
                        type1_error = reject_prob(n1, n2, r1, r, p0)
                        if type1_error > alpha:
                            continue

                        actual_power = reject_prob(n1, n2, r1, r, p1)
                        if actual_power < power:
                            continue

                        candidate = {
                            "n1": n1,
                            "n2": n2,
                            "n": n,
                            "r1": r1,
                            "r": r,
                            "pet0": pet0,
                            "en0": en0,
                            "actual_alpha": type1_error,
                            "actual_power": actual_power,
                        }

                        if (
                            best_opt is None
                            or candidate["en0"] < best_opt["en0"]
                            or (
                                abs(candidate["en0"] - best_opt["en0"]) < 1e-12
                                and candidate["n"] < best_opt["n"]
                            )
                        ):
                            best_opt = candidate

                        if (
                            best_min is None
                            or candidate["n"] < best_min["n"]
                            or (
                                candidate["n"] == best_min["n"]
                                and candidate["en0"] < best_min["en0"]
                            )
                        ):
                            best_min = candidate

                        # 对于固定(n,n1,r1)，r继续增大只会降低alpha和power，通常无需继续
                        break

        if best_opt is None or best_min is None:
            raise ValueError("在当前搜索范围内未找到满足约束的 Simon 二阶段设计。")

        selected = best_opt if design_type == "optimal" else best_min

        selected = dict(selected)
        selected["design_type"] = design_type
        selected["total_sample_size_with_dropout"] = int(
            math.ceil(selected["n"] / (1 - dropout_rate))
        )
        selected["n1_with_dropout"] = int(math.ceil(selected["n1"] / (1 - dropout_rate)))
        selected["n2_with_dropout"] = int(math.ceil(selected["n2"] / (1 - dropout_rate)))
        selected["optimal_design"] = {
            "n1": best_opt["n1"],
            "n2": best_opt["n2"],
            "n": best_opt["n"],
            "r1": best_opt["r1"],
            "r": best_opt["r"],
            "en0": best_opt["en0"],
            "pet0": best_opt["pet0"],
            "actual_alpha": best_opt["actual_alpha"],
            "actual_power": best_opt["actual_power"],
        }
        selected["minimax_design"] = {
            "n1": best_min["n1"],
            "n2": best_min["n2"],
            "n": best_min["n"],
            "r1": best_min["r1"],
            "r": best_min["r"],
            "en0": best_min["en0"],
            "pet0": best_min["pet0"],
            "actual_alpha": best_min["actual_alpha"],
            "actual_power": best_min["actual_power"],
        }
        return selected

    def _paired_nam_n(
        self,
        p10: float,
        p01: float,
        alpha: float,
        power: float,
        test_type: str,
        margin: float,
    ) -> tuple[float, float, float, float, float]:
        """Nam (1997) 型：返回 (n_raw, p_d, psi, eff, inner)。"""
        p_d = p10 + p01
        psi = p10 - p01
        inner = max(p_d - psi**2, 0.0)

        if test_type == "差异性检验":
            if abs(psi) < 1e-12:
                raise ValueError("π₁₀ 与 π₀₁ 相同，无法计算差异性样本量。")
            eff = abs(psi)
            z_a = float(norm.ppf(1 - alpha / 2))
            z_b = float(norm.ppf(power))
        elif test_type == "优效性检验":
            eff = psi - margin
            if eff <= 0:
                raise ValueError("优效性要求 π₁₀ − π₀₁ > Δ，请提高预期差异或缩小界值 Δ。")
            z_a = float(norm.ppf(1 - alpha))
            z_b = float(norm.ppf(power))
        elif test_type == "非劣效性检验":
            eff = psi + margin
            if eff <= 0:
                raise ValueError("非劣效性要求 π₁₀ − π₀₁ > −Δ，请检查预期与界值。")
            z_a = float(norm.ppf(1 - alpha))
            z_b = float(norm.ppf(power))
        elif test_type == "等效性检验":
            eff = margin - abs(psi)
            if eff <= 0:
                raise ValueError("等效性要求 |π₁₀ − π₀₁| < Δ，请调整预期或界值。")
            z_a = float(norm.ppf(1 - alpha))
            z_b = float(norm.ppf(1 - (1 - power) / 2.0))
        else:
            raise ValueError(f"不支持的检验类型: {test_type}")

        n_raw = (z_a * math.sqrt(p_d) + z_b * math.sqrt(inner)) ** 2 / (eff**2)
        return n_raw, p_d, psi, eff, inner

    def _paired_wald_blackwelder_n(
        self,
        p10: float,
        p01: float,
        alpha: float,
        power: float,
        test_type: str,
        margin: float,
    ) -> tuple[float, float]:
        """Blackwelder 型简化：忽略 (pd−ψ²) 修正，保守性略差，用于快速估算。"""
        p_d = p10 + p01
        psi = p10 - p01
        if test_type == "差异性检验":
            eff = abs(psi)
            z_sum = float(norm.ppf(1 - alpha / 2) + norm.ppf(power))
        elif test_type == "优效性检验":
            eff = psi - margin
            z_sum = float(norm.ppf(1 - alpha) + norm.ppf(power))
        elif test_type == "非劣效性检验":
            eff = psi + margin
            z_sum = float(norm.ppf(1 - alpha) + norm.ppf(power))
        elif test_type == "等效性检验":
            eff = margin - abs(psi)
            z_sum = float(norm.ppf(1 - alpha) + norm.ppf(1 - (1 - power) / 2.0))
        else:
            raise ValueError(f"不支持的检验类型: {test_type}")
        if eff <= 0:
            raise ValueError("当前效应量不足以满足该检验目标。")
        n_raw = (z_sum**2) * p_d / (eff**2)
        return n_raw, eff

    def _exact_mcnemar_diff_sample_size(
        self,
        p10: float,
        p01: float,
        alpha: float,
        power: float,
        dropout_rate: float,
    ) -> dict:
        """差异性：基于不一致对子集上二项精确检验（给定 n 时期望不一致对数）近似枚举。"""
        p_d = p10 + p01
        p_alt = p10 / p_d
        if p_d <= 0 or p_alt <= 0 or p_alt >= 1:
            raise ValueError("Exact McNemar 需要有效的 π₁₀、π₀₁。")

        def power_at_n(n: int) -> float:
            d = max(1, int(round(n * p_d)))
            reject: list[int] = []
            for k in range(d + 1):
                if float(binomtest(k, d, 0.5, alternative="two-sided").pvalue) <= alpha + 1e-15:
                    reject.append(k)
            return float(sum(binom.pmf(k, d, p_alt) for k in reject))

        n0 = int(max(5, math.ceil(self._paired_nam_n(p10, p01, alpha, power, "差异性检验", 0.0)[0])))
        n_start = max(3, n0 - 20)
        for n in range(n_start, 50001):
            if power_at_n(n) >= power - 1e-6:
                n_pairs = n
                n_drop = int(math.ceil(n_pairs / (1 - dropout_rate)))
                return {
                    "n_pairs": n_pairs,
                    "n_subjects": n_pairs,
                    "total_sample_size": n_pairs,
                    "n_pairs_with_dropout": n_drop,
                    "n_subjects_with_dropout": n_drop,
                    "total_sample_size_with_dropout": n_drop,
                    "method": "配对两组率：Exact McNemar（小样本 / discordant 很少，条件二项精确）",
                    "rate_note": (
                        "基于不一致对子集上双侧二项精确检验；与正态近似 McNemar 在大样本下趋于一致。"
                    ),
                    "paired_psi": float(p10 - p01),
                    "paired_pd": float(p_d),
                }
        raise ValueError("Exact McNemar 在可搜索范围内未找到满足效能的样本量。")

    def _exact_conditional_margin_sample_size(
        self,
        p10: float,
        p01: float,
        alpha: float,
        power: float,
        dropout_rate: float,
        test_type: str,
        margin: float,
    ) -> dict:
        """优效/非劣：不一致对子集上单侧二项精确检验，H₀ 取率差界值对应边界（补充/敏感性分析）。"""
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
            p0 = (p_d + margin) / (2.0 * p_d)
        else:
            if p_d <= margin + 1e-12:
                raise ValueError("需 pD > Δ，以保证 H₀ 边界有解。")
            p0 = (p_d - margin) / (2.0 * p_d)
        if not (0 < p0 < 1):
            raise ValueError("H₀ 边界下条件概率超出 (0,1)，请调整 pD 与 Δ。")

        def power_at_n(n: int) -> float:
            d = max(1, int(round(n * p_d)))
            reject: list[int] = []
            for k in range(d + 1):
                if float(binomtest(k, d, p0, alternative="greater").pvalue) <= alpha + 1e-15:
                    reject.append(k)
            return float(sum(binom.pmf(k, d, p_alt) for k in reject))

        n0 = int(
            max(
                5,
                math.ceil(
                    self._paired_nam_n(p10, p01, alpha, power, test_type, margin)[0]
                ),
            )
        )
        n_start = max(3, n0 - 40)
        for n in range(n_start, 50001):
            if power_at_n(n) >= power - 1e-6:
                n_drop = int(math.ceil(n / (1 - dropout_rate)))
                return {
                    "n_pairs": n,
                    "n_subjects": n,
                    "total_sample_size": n,
                    "n_pairs_with_dropout": n_drop,
                    "n_subjects_with_dropout": n_drop,
                    "total_sample_size_with_dropout": n_drop,
                    "n_treatment": n,
                    "n_control": 0,
                    "n_treatment_with_dropout": n_drop,
                    "n_control_with_dropout": 0,
                    "method": (
                        "配对两组率：exact conditional paired（不一致对条件二项单侧精确，H₀ 取界值边界）"
                    ),
                    "rate_note": (
                        "适用于极小样本或敏感性/补充分析；确证性主分析建议以 Nam Score 为主并交叉核对。"
                    ),
                    "paired_psi": float(p10 - p01),
                    "paired_pd": float(p_d),
                    "paired_margin": float(margin),
                }
        raise ValueError("exact conditional paired 在可搜索范围内未找到满足效能的样本量。")

    def compute_paired_proportion(
        self,
        p_discordant_treat: float,
        p_discordant_ctrl: float,
        alpha: float,
        power: float,
        dropout_rate: float,
        test_type: str,
        margin: float = 0.0,
        method: Optional[str] = None,
    ) -> dict:
        """
        配对两组率（McNemar / Nam Score / Wald Blackwelder / Exact）。

        π₁₀、π₀₁：不一致配对比例；margin：优效/非劣/等效的率差界值（与独立两组率同尺度）。

        method:
          - None / \"auto\"：差异性 → McNemar/Nam 若估计 n·p_d≥10 且 n≥25，否则 Exact McNemar；其余 → Nam Score
          - \"nam_score\"：Nam (1997) 型（差异性即 McNemar 正态近似）
          - \"wald_blackwelder\"：paired Wald-type 快速估算
          - \"exact_mcnemar\"：仅差异性
          - \"exact_conditional_paired\"：仅优效/非劣（条件二项精确）
        """
        p10, p01 = p_discordant_treat, p_discordant_ctrl
        if p10 <= 0 or p01 <= 0:
            raise ValueError("不一致配对比例（π₁₀、π₀₁）均须大于 0。")
        if p10 + p01 >= 1:
            raise ValueError("不一致配对比例之和须小于 1。")

        m = method or "auto"

        if m == "exact_conditional_paired":
            return self._exact_conditional_margin_sample_size(
                p10, p01, alpha, power, dropout_rate, test_type, margin
            )

        if m == "exact_mcnemar" and test_type != "差异性检验":
            raise ValueError("Exact McNemar 仅适用于差异性检验。")

        if (m == "auto" or m == "exact_mcnemar") and test_type == "差异性检验":
            if m == "exact_mcnemar":
                return self._exact_mcnemar_diff_sample_size(p10, p01, alpha, power, dropout_rate)
            n_raw, p_d, *_ = self._paired_nam_n(p10, p01, alpha, power, "差异性检验", 0.0)
            n_est = int(math.ceil(n_raw))
            if n_est * p_d >= 10.0 and n_est >= 25:
                m_use = "nam_score"
            else:
                return self._exact_mcnemar_diff_sample_size(p10, p01, alpha, power, dropout_rate)
        else:
            if m == "wald_blackwelder":
                m_use = "wald_blackwelder"
            else:
                m_use = "nam_score"

        if m_use == "wald_blackwelder":
            n_raw, eff = self._paired_wald_blackwelder_n(p10, p01, alpha, power, test_type, margin)
            inner_note = ""
            method_label = (
                "配对两组率：paired Wald-type approximation（Blackwelder 型简化，不建议作为主方法）"
            )
            rate_note = (
                "仅粗略估算；极端率或小样本下误差大，主分析请用 Nam Score 或 Exact 类方法复核。"
            )
        else:
            n_raw, p_d, psi, eff, inner = self._paired_nam_n(
                p10, p01, alpha, power, test_type, margin
            )
            inner_note = f"ψ={psi:.4f}，pD={p_d:.4f}。"
            if test_type == "等效性检验":
                method_label = "配对两组率：Nam score TOST（常规样本量，双单侧 Nam 型正态近似）"
                rate_note = (
                    "等效性（TOST）为双单侧检验；与独立样本 FM TOST 思路类似，基于 discordant pairs。"
                )
            elif test_type == "差异性检验":
                method_label = "配对两组率：McNemar 检验（正态近似，基于 discordant pairs；Nam 闭式）"
                rate_note = (
                    "配对二分类经典方法；小样本或 discordant 很少时请选 Exact McNemar。"
                )
            else:
                method_label = "配对两组率：Nam score 方法（优效/非劣，常规样本量）"
                rate_note = (
                    "配对二分类优/非劣常用标准近似；快速估算请用 Wald-type；极小样本可用 exact conditional paired 作敏感性分析。"
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
            "rate_note": rate_note,
            "paired_psi": float(p10 - p01),
            "paired_pd": float(p10 + p01),
            "paired_margin": float(margin),
        }
        if m_use == "nam_score" and test_type != "差异性检验":
            out["rate_note"] = rate_note + " " + inner_note
        elif m_use == "wald_blackwelder":
            out["rate_note"] = rate_note
        return out

    def compute_paired_mcnemar(
        self,
        p_discordant_treat: float,
        p_discordant_ctrl: float,
        alpha: float,
        power: float,
        dropout_rate: float,
        test_type: str,
    ) -> dict:
        """兼容旧接口：无界值 Nam 公式（差异性/未校正界值的优非劣等）。"""
        return self.compute_paired_proportion(
            p_discordant_treat,
            p_discordant_ctrl,
            alpha,
            power,
            dropout_rate,
            test_type,
            margin=0.0,
            method="nam_score",
        )

    @staticmethod
    def _three_arm_ps(p1: float, p2: float, p3: float) -> list:
        return [p1, p2, p3]

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

    def _three_arm_solve_n_ncx2(self, ncp_per_n: float, alpha: float, power: float, df: int) -> int:
        if ncp_per_n < 1e-12:
            raise ValueError("三组有效率相同或变换后效应量过小，无需计算样本量。")
        chi2_crit = chi2.ppf(1 - alpha, df)
        n = 5
        while n < 500000:
            pwr = 1 - ncx2.cdf(chi2_crit, df, n * ncp_per_n)
            if pwr >= power:
                return n
            n += 1
        raise ValueError("在合理范围内未能找到满足检验效能的样本量，请检查参数。")

    def _three_arm_mc_power(
        self,
        p1: float,
        p2: float,
        p3: float,
        n: int,
        alpha: float,
        n_reps: int,
        seed: int,
    ) -> float:
        """Pearson χ² 齐性检验（3×2）的蒙特卡洛功效（双侧，渐近 p 值）。"""
        rng = np.random.default_rng(seed)
        x1 = rng.binomial(n, p1, n_reps)
        x2 = rng.binomial(n, p2, n_reps)
        x3 = rng.binomial(n, p3, n_reps)
        hits = 0
        for i in range(n_reps):
            tab = np.array(
                [[float(x1[i]), float(n - x1[i])], [float(x2[i]), float(n - x2[i])], [float(x3[i]), float(n - x3[i])]],
                dtype=float,
            )
            try:
                _, p_val, _, _ = chi2_contingency(tab)
            except ValueError:
                continue
            if p_val < alpha:
                hits += 1
        return hits / float(n_reps)

    def _three_arm_mc_sample_size(
        self,
        p1: float,
        p2: float,
        p3: float,
        alpha: float,
        power: float,
        n_reps: int,
        seed: int,
        n_hint: int,
    ) -> tuple:
        """搜索满足蒙特卡洛功效的每组最小 n；从提示值附近开始以控制耗时。"""
        n = max(5, n_hint - 30)
        last_pw = 0.0
        while n < 500000:
            last_pw = self._three_arm_mc_power(p1, p2, p3, n, alpha, n_reps, seed + n % 997)
            if last_pw >= power - 0.008:
                break
            if n < 80:
                n += 1
            elif n < 500:
                n += 2
            else:
                n += max(5, int(n * 0.04))
        else:
            raise ValueError("蒙特卡洛搜索在可及范围内未找到满足效能的样本量，请放宽效应或改用渐近法。")
        while n > 5:
            pw_prev = self._three_arm_mc_power(p1, p2, p3, n - 1, alpha, n_reps, seed + (n - 1) % 997)
            if pw_prev < power - 0.008:
                break
            n -= 1
        return n, last_pw

    @staticmethod
    def _cochran_armitage_pvalue_two_sided(x: np.ndarray, n_per_group: int, scores: np.ndarray) -> float:
        """Cochran–Armitage 趋势检验，双侧渐近 p 值（等额分组）。"""
        x = np.asarray(x, dtype=float)
        scores = np.asarray(scores, dtype=float)
        k = len(x)
        n_i = np.full(k, float(n_per_group), dtype=float)
        n_total = float(np.sum(n_i))
        if n_total <= 0:
            return 1.0
        p_hat = float(np.sum(x) / n_total)
        t_bar = float(np.sum(n_i * scores) / n_total)
        t_stat = float(np.sum((scores - t_bar) * (x - n_i * p_hat)))
        var_t = p_hat * (1.0 - p_hat) * float(np.sum(n_i * (scores - t_bar) ** 2))
        if var_t <= 1e-15:
            return 1.0
        z = t_stat / math.sqrt(var_t)
        p = 2.0 * (1.0 - float(norm.cdf(abs(z))))
        return float(min(max(p, 0.0), 1.0))

    def _three_arm_ca_mc_power(
        self,
        p1: float,
        p2: float,
        p3: float,
        n: int,
        alpha: float,
        scores: np.ndarray,
        n_reps: int,
        seed: int,
    ) -> float:
        rng = np.random.default_rng(seed)
        x1 = rng.binomial(n, p1, n_reps)
        x2 = rng.binomial(n, p2, n_reps)
        x3 = rng.binomial(n, p3, n_reps)
        hits = 0
        for i in range(n_reps):
            x = np.array([float(x1[i]), float(x2[i]), float(x3[i])], dtype=float)
            p_val = self._cochran_armitage_pvalue_two_sided(x, n, scores)
            if p_val < alpha:
                hits += 1
        return hits / float(n_reps)

    def _three_arm_ca_mc_sample_size(
        self,
        p1: float,
        p2: float,
        p3: float,
        alpha: float,
        power: float,
        scores: np.ndarray,
        n_reps: int,
        seed: int,
        n_hint: int,
    ) -> tuple:
        n = max(5, n_hint - 20)
        last_pw = 0.0
        while n < 500000:
            last_pw = self._three_arm_ca_mc_power(p1, p2, p3, n, alpha, scores, n_reps, seed + n % 997)
            if last_pw >= power - 0.008:
                break
            if n < 80:
                n += 1
            elif n < 500:
                n += 2
            else:
                n += max(5, int(n * 0.04))
        else:
            raise ValueError("蒙特卡洛搜索在可及范围内未找到满足效能的样本量，请放宽效应或检查趋势假设。")
        while n > 5:
            pw_prev = self._three_arm_ca_mc_power(p1, p2, p3, n - 1, alpha, scores, n_reps, seed + (n - 1) % 997)
            if pw_prev < power - 0.008:
                break
            n -= 1
        return n, last_pw

    def compute_three_arm_cochran_armitage(
        self,
        p1: float,
        p2: float,
        p3: float,
        alpha: float,
        power: float,
        dropout_rate: float,
        scores: tuple = (0.0, 1.0, 2.0),
        mc_replications: int = 8000,
        mc_seed: int = 42,
    ) -> dict:
        """三组率 Cochran–Armitage 线性趋势检验（双侧），样本量由蒙特卡洛估计（等额分配）。"""
        ps = self._three_arm_ps(p1, p2, p3)
        k = 3
        sc = np.array(scores, dtype=float)
        if len(sc) != k:
            raise ValueError("趋势分数须为 3 个实数，与三组一一对应。")
        if not (sc[0] < sc[1] < sc[2]):
            raise ValueError("趋势分数须严格递增（例如 0→1→2，对应剂量/序次）。")
        p_bar = sum(ps) / k
        if p_bar <= 1e-9 or p_bar >= 1 - 1e-9:
            raise ValueError("三组平均有效率过于极端，无法计算。")
        n_i = np.full(k, 1.0, dtype=float)
        n_total = float(k)
        t_bar = float(np.sum(n_i * sc) / n_total)
        trend_contrast = sum((sc[i] - t_bar) * (ps[i] - p_bar) for i in range(k))
        if abs(trend_contrast) < 1e-12:
            raise ValueError("在假定分数下，三组率相对均值的线性趋势分量为 0，无需按趋势检验计算样本量。")

        ncp_p = self._three_arm_ncp_pearson(ps)
        if ncp_p >= 1e-12:
            n_pearson = self._three_arm_solve_n_ncx2(ncp_p, alpha, power, k - 1)
            n_hint = min(max(5, n_pearson), 3000)
        else:
            n_hint = 50
        n, mc_pw = self._three_arm_ca_mc_sample_size(
            p1, p2, p3, alpha, power, sc, mc_replications, mc_seed, n_hint
        )
        n_per_group_drop = int(math.ceil(n / (1 - dropout_rate)))
        min_cell = self._three_arm_min_expected_cell(ps, n)
        return {
            "n_per_group": n,
            "total_sample_size": k * n,
            "n_per_group_with_dropout": n_per_group_drop,
            "total_sample_size_with_dropout": k * n_per_group_drop,
            "method": "三组独立样本率：Cochran–Armitage 趋势检验（蒙特卡洛双侧功效，等额分配）",
            "rate_note": (
                f"假定三组有序且关注线性趋势；分数为 ({scores[0]:g}, {scores[1]:g}, {scores[2]:g})。"
                "以渐近 Z 统计量双侧 p 值做蒙特卡洛功效；大样本可与正态近似对照。注册试验请统计师复核。"
            ),
            "three_arm_method": "cochran_armitage_mc",
            "min_expected_cell_approx": float(min_cell),
            "actual_power_mc": float(mc_pw),
            "mc_replications": int(mc_replications),
            "ca_scores": (float(sc[0]), float(sc[1]), float(sc[2])),
        }

    def compute_three_arm_proportion(
        self,
        p1: float,
        p2: float,
        p3: float,
        alpha: float,
        power: float,
        dropout_rate: float,
        method: Optional[str] = None,
        mc_replications: int = 8000,
        mc_seed: int = 42,
    ) -> dict:
        """三组率差异性检验样本量（k×2 齐性，等额分配，df=2）。

        method: None=自动（期望频数充分→Pearson，否则→蒙特卡洛）；pearson；monte_carlo。
        """
        ps = self._three_arm_ps(p1, p2, p3)
        k = 3
        df = k - 1

        p_bar = sum(ps) / k
        if p_bar <= 1e-9 or p_bar >= 1 - 1e-9:
            raise ValueError("三组平均有效率过于极端，无法计算。")

        ncp_p = self._three_arm_ncp_pearson(ps)
        if ncp_p < 1e-12:
            raise ValueError("三组有效率相同或效应量过小，无需计算样本量。")

        m = method or "auto"
        if m == "auto":
            normal_zone = all(0.2 <= pi <= 0.8 for pi in ps)
            n_pearson = self._three_arm_solve_n_ncx2(ncp_p, alpha, power, df)
            min_cell_p = self._three_arm_min_expected_cell(ps, n_pearson)
            if normal_zone and min_cell_p >= 5.0 - 1e-9:
                m = "pearson"
            else:
                m = "monte_carlo"

        if m == "pearson":
            ncp_use = ncp_p
            n = self._three_arm_solve_n_ncx2(ncp_use, alpha, power, df)
            method_label = "三组独立样本率：Pearson χ² 齐性检验（非中心 χ² 法，等额分配）"
            rate_note = (
                "基于渐近 Pearson χ² 的非中心参数；与合并方差思路一致，适用于各组期望频数较充分时。"
                "若率极偏或每组 n 很小，请查看「自动」或手动选蒙特卡洛。"
            )
        elif m == "monte_carlo":
            n_hint = self._three_arm_solve_n_ncx2(ncp_p, alpha, power, df)
            n, mc_pw = self._three_arm_mc_sample_size(
                p1, p2, p3, alpha, power, mc_replications, mc_seed, n_hint
            )
            method_label = "三组独立样本率：蒙特卡洛模拟（Pearson χ² 检验统计量，小样本/复核）"
            rate_note = (
                f"在每组 n 下用二项随机数模拟 3×2 表，以 Pearson χ² 渐近 p 值估计功效（重复 {mc_replications} 次）。"
                "注册试验建议与统计师复核；大样本时与非中心 χ² 接近。"
            )
            ncp_use = float("nan")
        else:
            raise ValueError(f"不支持的三组率方法: {m}")

        n_per_group_drop = int(math.ceil(n / (1 - dropout_rate)))
        min_cell = self._three_arm_min_expected_cell(ps, n)
        out: Dict[str, Any] = {
            "n_per_group": n,
            "total_sample_size": k * n,
            "n_per_group_with_dropout": n_per_group_drop,
            "total_sample_size_with_dropout": k * n_per_group_drop,
            "method": method_label,
            "rate_note": rate_note,
            "three_arm_method": m,
            "min_expected_cell_approx": float(min_cell),
        }
        if m != "monte_carlo":
            out["ncp_per_n"] = float(ncp_use)
        else:
            out["actual_power_mc"] = float(mc_pw)
            out["mc_replications"] = int(mc_replications)
        return out

    def _format_result(self, nt: float, nc: float, dropout: float) -> dict:
        nt_ceil, nc_ceil = int(math.ceil(nt)), int(math.ceil(nc))
        nt_drop, nc_drop = int(math.ceil(nt_ceil / (1 - dropout))), int(math.ceil(nc_ceil / (1 - dropout)))
        return {
            "n_treatment": nt_ceil, "n_control": nc_ceil, "total_sample_size": nt_ceil + nc_ceil,
            "n_treatment_with_dropout": nt_drop, "n_control_with_dropout": nc_drop, "total_sample_size_with_dropout": nt_drop + nc_drop
        }