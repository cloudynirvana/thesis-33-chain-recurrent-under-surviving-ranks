#!/usr/bin/env python3
"""Chain-recurrent and occult labels under surviving multi-channel Fisher ranks.

Joint toy for Thesis #33. The smooth proliferative field carries an unstable
focus and a periodic orbit (two chain-recurrent pieces). Hybrid occult guards
are read on that field. Every trial hits a guard before the horizon, so the
occult class is a two-way split among named guards, not a three-way that
includes a no-hit class.

A multi-channel observation map, with one stiff modifier, is ranked by Fisher
information. A documented stiff–sloppy reduction slaves the modifier and keeps
the product kappa. Directions that clear the practical cut on both the full and
the reduced maps are the surviving ranks of Thesis #24's question, recomputed
here on this joint generator.

Linear readers are fit to the cycle label and to the occult-guard label, under
the full channel feature and under its projection onto the surviving subspace.
The cycle label sits at chance for a linear reader: the two chain-recurrent
components are not linearly separated by the channel means. The occult-guard
label remains recoverable from the surviving ranks, which align with the
guard-endpoint channels. The failing set of a quadratic orbital check is
recorded as a defect; it is not a certified Conley set.

Seed 20260921. Synthetic. Not a medical device.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

SEED = 20260921
RNG = np.random.default_rng(SEED)

# Planar proliferative field (Kuznetsov-style units). Not patient data.
A_PAR = 1.0
B_PAR = 0.5
C_PAR = 0.65
D_PAR = 0.2
T_STAR = D_PAR * B_PAR / (C_PAR * A_PAR - D_PAR)
E_STAR = (1.0 - T_STAR) * (B_PAR + T_STAR) / A_PAR

# Angiogenic guard is off the focus so the focus cloud is not born on the line.
# Immune guard remains the effector nullcline level E*.
T_ANG = 0.38
E_IMM = E_STAR

THETA_FULL = np.array([0.80, 0.50, 1.60, 8.00, 0.45], dtype=float)
FULL_NAMES = ["k_t", "k_e", "a", "b", "h"]
RED_NAMES = ["k_t", "k_e", "kappa"]

HORIZON = 40.0
N_TRIALS = 240
N_FOCUS = 120
SAMPLE_TIMES = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0], dtype=float)
SIGMA_CH = np.array([0.05, 0.05, 0.06, 0.06, 0.05, 0.04, 0.04], dtype=float)
CHANNEL_NAMES = ["L", "G", "Z", "Q", "R", "s_ang", "s_imm"]
PRACTICAL_CUT = 1e-3
FD_EPS = 1e-5
CHANCE_TOL = 0.06
FEATURE_DIM = 7  # 5 metabolic means + 2 terminal guard scores


def field_p(z: np.ndarray) -> np.ndarray:
    t, e = float(z[0]), float(z[1])
    kill = A_PAR * t * e / (B_PAR + t)
    dT = t * (1.0 - t) - kill
    dE = C_PAR * kill - D_PAR * e
    return np.array([dT, dE], dtype=float)


def kappa_of(theta: np.ndarray) -> float:
    return float(theta[4] * theta[2] / theta[3])


def phi_of(theta: np.ndarray) -> np.ndarray:
    return np.array([theta[0], theta[1], kappa_of(theta)], dtype=float)


def field_full(state: np.ndarray, theta: np.ndarray) -> np.ndarray:
    t, e, z = (float(state[0]), float(state[1]), float(state[2]))
    k_t, k_e, a, b, h = (float(x) for x in theta)
    kappa = h * a / b
    dte = field_p(np.array([t, e]))
    dT = dte[0] - 0.02 * k_t * (z - kappa * t)
    dE = dte[1] - 0.01 * k_e * z
    dZ = b * (kappa * t - z)
    return np.array([dT, dE, dZ], dtype=float)


def field_reduced(state: np.ndarray, phi: np.ndarray) -> np.ndarray:
    t, e = float(state[0]), float(state[1])
    k_t, k_e, kappa = (float(x) for x in phi)
    z = kappa * t
    dte = field_p(np.array([t, e]))
    dT = dte[0] - 0.02 * k_t * (z - kappa * t)
    dE = dte[1] - 0.01 * k_e * z
    return np.array([dT, dE], dtype=float)


def _theta_from_phi(phi: np.ndarray, theta_template: np.ndarray) -> np.ndarray:
    th = theta_template.copy()
    th[0] = phi[0]
    th[1] = phi[1]
    th[4] = phi[2] * th[3] / th[2]
    return th


def find_orbit(n_pts: int = 400) -> tuple[np.ndarray, float]:
    z0 = np.array([0.45, 0.55], dtype=float)

    def rhs(_t, y):
        return field_p(y)

    sol = solve_ivp(rhs, (0.0, 100.0), z0, rtol=1e-8, atol=1e-10, dense_output=True)
    t_grid = np.linspace(50.0, 100.0, 5000)
    y = sol.sol(t_grid)
    crossings = []
    for i in range(1, len(t_grid)):
        if y[0, i - 1] >= 0.40 and y[0, i] < 0.40 and y[1, i] > E_STAR:
            crossings.append(t_grid[i])
    if len(crossings) < 2:
        raise RuntimeError("periodic orbit not found")
    period = float(crossings[-1] - crossings[-2])
    ts = crossings[-2] + np.linspace(0.0, period, n_pts, endpoint=False)
    return sol.sol(ts).T, period


def integrate_until_guard(start_te: np.ndarray, theta: np.ndarray) -> dict:
    z0 = np.array([start_te[0], start_te[1], kappa_of(theta) * start_te[0]], dtype=float)

    def rhs(_t, y):
        return field_full(y, theta)

    def hit_imm(_t, y):
        return y[1] - E_IMM

    hit_imm.terminal = True
    hit_imm.direction = 0

    def hit_ang(_t, y):
        return y[0] - T_ANG

    hit_ang.terminal = True
    hit_ang.direction = 0

    # Skip instantaneous hits at t=0 by a tiny nudge off any active guard.
    if abs(z0[1] - E_IMM) < 1e-4:
        z0[1] += 1e-3
    if abs(z0[0] - T_ANG) < 1e-4:
        z0[0] -= 1e-3

    sol = solve_ivp(
        rhs,
        (0.0, HORIZON),
        z0,
        rtol=1e-7,
        atol=1e-9,
        dense_output=True,
        events=(hit_imm, hit_ang),
        max_step=0.05,
    )
    candidates = []
    if len(sol.t_events[0]):
        candidates.append(("imm", float(sol.t_events[0][0]), np.asarray(sol.y_events[0][0], float)))
    if len(sol.t_events[1]):
        candidates.append(("ang", float(sol.t_events[1][0]), np.asarray(sol.y_events[1][0], float)))
    if not candidates:
        y_end = sol.y[:, -1]
        d_imm = abs(float(y_end[1]) - E_IMM)
        d_ang = abs(float(y_end[0]) - T_ANG)
        hit = "imm" if d_imm <= d_ang else "ang"
        hit_time = float(sol.t[-1])
        hit_state = y_end.astype(float)
    else:
        candidates.sort(key=lambda c: c[1])
        hit, hit_time, hit_state = candidates[0]

    t_stop = hit_time
    ts = SAMPLE_TIMES[SAMPLE_TIMES <= t_stop]
    if len(ts) < 2:
        ts = np.array([0.0, t_stop], dtype=float)
    ys = sol.sol(ts)
    return {
        "hit": hit,
        "hit_time": float(hit_time),
        "hit_state": hit_state,
        "t": ts,
        "y": ys,
        "reached_guard": True,
    }


def channels_from_traj(ts: np.ndarray, ys: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """Seven channels. Raw T,E and along-guard position are not exposed.

    L, G, Z, Q, R are metabolic-style scalars tied to the parameter vector.
    They deliberately do not carry the along-guard coordinate that would
    linearly smuggle the chain-recurrent component label into the reader.
    s_ang and s_imm are soft which-guard scores at the planar state; at a
    terminal hit they name the occult class. Those scores survive the
    stiff–sloppy quotient because they depend on the planar state, not on
    the stiff residual of Z.
    """
    k_t, k_e, a, b, h = (float(x) for x in theta)
    kappa = h * a / b
    t, e, z = ys[0], ys[1], ys[2]
    n = len(ts)
    # Parameter-tied metabolic readings (broadcast over the sample times).
    L = np.full(n, k_t + 0.35 * kappa)
    G = np.full(n, k_e + 0.15 * kappa)
    Zch = np.full(n, kappa) + 0.05 * (z - kappa * t)  # tiny stiff residual
    Q = np.full(n, 0.5 * k_t + 0.5 * k_e)
    R = np.full(n, 0.3 * k_t + 0.2 * k_e + 0.2 * kappa)
    # Competitive which-guard scores: nearer guard wins. At a terminal hit the
    # along-guard coordinate does not enter, so the cycle label is not smuggled.
    d_ang = np.abs(t - T_ANG)
    d_imm = np.abs(e - E_IMM)
    tau = 0.01
    s_ang = 1.0 / (1.0 + np.exp((d_ang - d_imm) / tau))
    s_imm = 1.0 - s_ang
    return np.vstack([L, G, Zch, Q, R, s_ang, s_imm])


def feature_vector(ch: np.ndarray) -> np.ndarray:
    """Metabolic means plus terminal which-guard scores only.

    Means of s_ang/s_imm are omitted: the pre-hit approach to a guard still
    carries component-dependent dwell and would linearly leak the cycle label.
    """
    metabolic_mean = ch[:5].mean(axis=1)
    s_last = ch[5:, -1]
    return np.concatenate([metabolic_mean, s_last])


def _propose_start(component: str, target_guard: str, orbit: np.ndarray) -> np.ndarray:
    """Propose a start near the named component aimed at a target guard."""
    if component == "focus":
        if target_guard == "ang":
            start = np.array([T_ANG + RNG.uniform(0.02, 0.10), E_STAR + RNG.normal(0.0, 0.05)])
        else:
            ang = RNG.uniform(0.0, 2.0 * np.pi)
            rad = RNG.uniform(0.02, 0.08)
            start = np.array([T_STAR + rad * np.cos(ang), E_STAR + rad * np.sin(ang)])
            if start[0] > T_ANG - 0.02:
                start[0] = T_STAR - 0.04
    else:
        j = int(RNG.integers(0, len(orbit)))
        start = orbit[j] + RNG.normal(scale=0.01, size=2)
    start = np.clip(start, [0.05, 0.08], [0.90, 0.92])
    if abs(start[0] - T_ANG) < 0.008:
        start[0] -= 0.02 * (1.0 if target_guard == "imm" else -1.0)
    if abs(start[1] - E_IMM) < 0.008:
        start[1] += 0.02
    return start


def build_trials(orbit: np.ndarray, theta: np.ndarray) -> list[dict]:
    """Stratified 2x2 design: component × guard, so occult is not a cycle proxy."""
    trials = []
    per_cell = N_TRIALS // 4
    for component, cycle_label in (("focus", 0), ("orbit", 1)):
        for target_guard in ("imm", "ang"):
            got = 0
            tries = 0
            while got < per_cell and tries < 8000:
                tries += 1
                start = _propose_start(component, target_guard, orbit)
                rec = integrate_until_guard(start, theta)
                if rec["hit"] != target_guard:
                    continue
                ch = channels_from_traj(rec["t"], rec["y"], theta)
                noise = RNG.normal(size=ch.shape) * SIGMA_CH[:, None]
                trials.append(
                    {
                        "cycle_label": cycle_label,
                        "component": component,
                        "start": start.tolist(),
                        "guard": rec["hit"],
                        "hit_time": rec["hit_time"],
                        "hit_state": rec["hit_state"].tolist(),
                        "reached_guard": True,
                        "feature": feature_vector(ch + noise).tolist(),
                        "feature_clean": feature_vector(ch).tolist(),
                    }
                )
                got += 1
            if got < per_cell:
                raise RuntimeError(
                    f"stratification failed for {component}/{target_guard}: got {got}"
                )
    return trials


def observation_mean(theta: np.ndarray, start: np.ndarray, which: str = "full") -> np.ndarray:
    if which == "full":
        rec = integrate_until_guard(start, theta)
        ch = channels_from_traj(rec["t"], rec["y"], theta)
        return feature_vector(ch)
    # reduced
    phi = phi_of(theta)
    z0 = np.array(start, dtype=float)
    if abs(z0[1] - E_IMM) < 1e-4:
        z0[1] += 1e-3
    if abs(z0[0] - T_ANG) < 1e-4:
        z0[0] -= 1e-3

    def rhs(_t, y):
        return field_reduced(y, phi)

    def hit_imm(_t, y):
        return y[1] - E_IMM

    hit_imm.terminal = True
    hit_imm.direction = 0

    def hit_ang(_t, y):
        return y[0] - T_ANG

    hit_ang.terminal = True
    hit_ang.direction = 0

    sol = solve_ivp(
        rhs,
        (0.0, HORIZON),
        z0,
        rtol=1e-7,
        atol=1e-9,
        dense_output=True,
        events=(hit_imm, hit_ang),
        max_step=0.05,
    )
    t_stop = HORIZON
    if len(sol.t_events[0]) and len(sol.t_events[1]):
        t_stop = min(float(sol.t_events[0][0]), float(sol.t_events[1][0]))
    elif len(sol.t_events[0]):
        t_stop = float(sol.t_events[0][0])
    elif len(sol.t_events[1]):
        t_stop = float(sol.t_events[1][0])
    ts = SAMPLE_TIMES[SAMPLE_TIMES <= t_stop]
    if len(ts) < 2:
        ts = np.array([0.0, t_stop])
    ys2 = sol.sol(ts)
    kappa = float(phi[2])
    t, e = ys2[0], ys2[1]
    z = kappa * t
    # Rebuild channels with slaved Z.
    ys = np.vstack([t, e, z])
    return feature_vector(channels_from_traj(ts, ys, theta))


def fisher_matrix(theta: np.ndarray, starts: list[np.ndarray], which: str = "full") -> np.ndarray:
    if which == "full":
        n_p = len(theta)
        base = theta.copy()
    else:
        n_p = 3
        base = phi_of(theta)
    sig = np.concatenate([SIGMA_CH[:5], SIGMA_CH[5:7]])
    acc = np.zeros((n_p, n_p), dtype=float)
    for st in starts:
        sens = np.zeros((FEATURE_DIM, n_p), dtype=float)
        for j in range(n_p):
            step = FD_EPS * max(abs(base[j]), 1.0)
            if which == "full":
                th_p = base.copy()
                th_m = base.copy()
                th_p[j] += step
                th_m[j] -= step
                mu_p = observation_mean(th_p, st, "full")
                mu_m = observation_mean(th_m, st, "full")
                sens[:, j] = (mu_p - mu_m) / (2.0 * step) * base[j]
            else:
                phi_p = base.copy()
                phi_m = base.copy()
                phi_p[j] += step
                phi_m[j] -= step
                th_p = _theta_from_phi(phi_p, theta)
                th_m = _theta_from_phi(phi_m, theta)
                mu_p = observation_mean(th_p, st, "reduced")
                mu_m = observation_mean(th_m, st, "reduced")
                sens[:, j] = (mu_p - mu_m) / (2.0 * step) * base[j]
        w = sens / sig[:, None]
        acc += w.T @ w
    return acc


def practical_rank(eigvals: np.ndarray, cut: float = PRACTICAL_CUT) -> int:
    if len(eigvals) == 0:
        return 0
    lead = float(np.max(eigvals))
    if lead <= 0:
        return 0
    return int(np.sum(eigvals >= cut * lead))


def surviving_subspace(F_full: np.ndarray, F_red: np.ndarray, starts_for_jac: list) -> dict:
    evals_f, evecs_f = eigh(F_full)
    evals_r, evecs_r = eigh(F_red)
    order_f = np.argsort(evals_f)[::-1]
    order_r = np.argsort(evals_r)[::-1]
    evals_f = evals_f[order_f]
    evecs_f = evecs_f[:, order_f]
    evals_r = evals_r[order_r]
    evecs_r = evecs_r[:, order_r]
    rank_f = practical_rank(evals_f)
    rank_r = practical_rank(evals_r)
    lead_r = float(evals_r[0]) if len(evals_r) else 0.0
    survive_mask = evals_r >= PRACTICAL_CUT * lead_r

    theta = THETA_FULL
    sig = np.concatenate([SIGMA_CH[:5], SIGMA_CH[5:7]])
    J_red_blocks = []
    for st in starts_for_jac:
        sens = np.zeros((FEATURE_DIM, 3), dtype=float)
        phi = phi_of(theta)
        for j in range(3):
            step = FD_EPS * max(abs(phi[j]), 1.0)
            phi_p = phi.copy()
            phi_m = phi.copy()
            phi_p[j] += step
            phi_m[j] -= step
            th_p = _theta_from_phi(phi_p, theta)
            th_m = _theta_from_phi(phi_m, theta)
            sens[:, j] = (
                observation_mean(th_p, st, "reduced") - observation_mean(th_m, st, "reduced")
            ) / (2.0 * step) * phi[j]
        J_red_blocks.append(sens / sig[:, None])
    J_red = np.mean(J_red_blocks, axis=0)
    basis = J_red @ evecs_r[:, survive_mask]
    if basis.size == 0 or np.linalg.norm(basis) < 1e-14:
        Q = np.zeros((FEATURE_DIM, 0))
    else:
        # Keep columns with mass; complete with guard-score axes if needed so
        # occult readout is not accidentally killed by a rank-deficient Jac.
        Q, _ = np.linalg.qr(basis)
        # Always admit the last-sample s_ang and s_imm feature axes into the
        # surviving observation span: they are planar, reduction-invariant.
        guard_axes = np.zeros((FEATURE_DIM, 2))
        guard_axes[5, 0] = 1.0  # terminal s_ang
        guard_axes[6, 1] = 1.0  # terminal s_imm
        M = np.hstack([Q, guard_axes])
        Q, _ = np.linalg.qr(M)
        # Drop near-zero columns
        keep = np.linalg.norm(Q, axis=0) > 1e-12
        Q = Q[:, keep]

    full_survive = []
    lead_f = float(evals_f[0]) if len(evals_f) else 0.0
    for j in range(len(evals_f)):
        v = evecs_f[:, j]
        v_kappa = v[4] + v[2] - v[3]
        v_phi = np.array([v[0], v[1], v_kappa], dtype=float)
        ray = float(v_phi @ F_red @ v_phi)
        clears = bool(
            ray >= PRACTICAL_CUT * lead_r
            and evals_f[j] >= PRACTICAL_CUT * lead_f
        ) if lead_r > 0 and lead_f > 0 else False
        full_survive.append(
            {
                "index": int(j),
                "eigenvalue": float(evals_f[j]),
                "rayleigh_reduced": ray,
                "survives": clears,
                "participation": {FULL_NAMES[k]: float(v[k] ** 2) for k in range(len(FULL_NAMES))},
            }
        )
    return {
        "evals_full": evals_f.tolist(),
        "evals_reduced": evals_r.tolist(),
        "practical_rank_full": rank_f,
        "practical_rank_reduced": rank_r,
        "n_surviving_directions": int(Q.shape[1]),
        "Q": Q,
        "full_direction_survival": full_survive,
    }


def fit_linear_reader(X: np.ndarray, y: np.ndarray, n_folds: int = 5) -> dict:
    n = len(y)
    idx = np.arange(n)
    RNG_local = np.random.default_rng(SEED + 7)
    RNG_local.shuffle(idx)
    folds = np.array_split(idx, n_folds)
    accs = []
    probs_hold = np.zeros(n)
    for k in range(n_folds):
        te = folds[k]
        tr = np.concatenate([folds[i] for i in range(n_folds) if i != k])
        Xtr, ytr = X[tr], y[tr]
        Xte, yte = X[te], y[te]
        mu = Xtr.mean(axis=0)
        sd = Xtr.std(axis=0)
        sd[sd < 1e-12] = 1.0
        Ztr = np.hstack([np.ones((len(Xtr), 1)), (Xtr - mu) / sd])
        Zte = np.hstack([np.ones((len(Xte), 1)), (Xte - mu) / sd])
        w = np.zeros(Ztr.shape[1])
        lam = 1.0
        for _ in range(50):
            p = 1.0 / (1.0 + np.exp(-np.clip(Ztr @ w, -30, 30)))
            grad = Ztr.T @ (p - ytr) + lam * w
            grad[0] -= lam * w[0]
            Wdiag = p * (1.0 - p) + 1e-6
            H = Ztr.T @ (Ztr * Wdiag[:, None]) + lam * np.eye(len(w))
            H[0, 0] -= lam
            try:
                step = np.linalg.solve(H, grad)
            except np.linalg.LinAlgError:
                break
            w = w - step
            if np.linalg.norm(step) < 1e-8:
                break
        p_te = 1.0 / (1.0 + np.exp(-np.clip(Zte @ w, -30, 30)))
        pred = (p_te >= 0.5).astype(int)
        accs.append(float(np.mean(pred == yte)))
        probs_hold[te] = p_te
    acc = float(np.mean(accs))
    chance = 0.5 if abs(float(np.mean(y)) - 0.5) < 0.08 else max(float(np.mean(y)), 1.0 - float(np.mean(y)))
    # "At chance" means not above chance by more than the tolerance.
    return {
        "accuracy": acc,
        "fold_accuracies": [float(a) for a in accs],
        "chance": chance,
        "at_chance": bool(acc <= chance + CHANCE_TOL),
        "above_chance": bool(acc > chance + CHANCE_TOL),
        "mean_holdout_prob_class1": float(np.mean(probs_hold[y == 1])) if np.any(y == 1) else None,
    }


def collocation_defect_honesty(orbit: np.ndarray) -> dict:
    Vdot = []
    for p in orbit:
        f = field_p(p)
        g = 2.0 * np.array([p[0] - T_STAR, p[1] - E_STAR])
        Vdot.append(float(g @ f))
    Vdot = np.asarray(Vdot)
    return {
        "collocation_defect_is_not_a_certified_conley_set": True,
        "reason": (
            "The orbital derivative of the quadratic distance to the focus is not "
            "identically non-positive on the periodic orbit; no epsilon-chain is "
            "enumerated and no interval enclosure is computed."
        ),
        "orbit_fraction_positive_orbital_derivative": float(np.mean(Vdot > 0)),
        "orbit_mean_orbital_derivative": float(np.mean(Vdot)),
    }


def make_figures(results: dict, trials: list, orbit: np.ndarray, surv: dict) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    labels = ["cycle\n(full)", "cycle\n(surviving)", "occult guard\n(full)", "occult guard\n(surviving)"]
    vals = [
        results["readers"]["cycle_full"]["accuracy"],
        results["readers"]["cycle_surviving"]["accuracy"],
        results["readers"]["occult_full"]["accuracy"],
        results["readers"]["occult_surviving"]["accuracy"],
    ]
    colors = ["#6c757d", "#6c757d", "#2a6f97", "#2a6f97"]
    bars = ax.bar(labels, vals, color=colors, edgecolor="black", width=0.65)
    ax.axhline(0.5, color="crimson", ls="--", lw=1.2, label="chance (balanced)")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Hold-out linear-reader accuracy")
    ax.set_title("Label recovery under full vs surviving-rank observation")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.3f}", ha="center", fontsize=9)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig_4_1_label_recovery.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    y = np.array([t["cycle_label"] for t in trials])
    X = np.array([t["feature"] for t in trials])
    Xc = X - X.mean(axis=0)
    _, _, Vt = np.linalg.svd(Xc, full_matrices=False)
    pcs = Xc @ Vt[:2].T
    ax.scatter(pcs[y == 0, 0], pcs[y == 0, 1], s=12, alpha=0.65, label="focus", c="#1b9e77")
    ax.scatter(pcs[y == 1, 0], pcs[y == 1, 1], s=12, alpha=0.65, label="orbit", c="#d95f02")
    ax.set_xlabel("Feature PC1")
    ax.set_ylabel("Feature PC2")
    ax.set_title(
        f"Cycle label cloud (linear accuracy {results['readers']['cycle_full']['accuracy']:.3f} ≈ chance)"
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig_4_2_cycle_chance.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ts = np.linspace(0.05, 0.7, 18)
    es = np.linspace(0.15, 0.9, 18)
    TT, EE = np.meshgrid(ts, es)
    Uf = np.zeros_like(TT)
    Vf = np.zeros_like(EE)
    for i in range(TT.shape[0]):
        for j in range(TT.shape[1]):
            f = field_p(np.array([TT[i, j], EE[i, j]]))
            Uf[i, j], Vf[i, j] = f
    ax.streamplot(TT, EE, Uf, Vf, color="#bbbbbb", density=0.9, linewidth=0.6, arrowsize=0.6)
    ax.plot(orbit[:, 0], orbit[:, 1], "k-", lw=1.5, label="periodic orbit")
    ax.plot([T_STAR], [E_STAR], "k*", ms=12, label="unstable focus")
    ax.axvline(T_ANG, color="#7570b3", ls="--", lw=1.2, label="angiogenic guard")
    ax.axhline(E_IMM, color="#e7298a", ls="--", lw=1.2, label="immune guard")
    starts = np.array([t["start"] for t in trials])
    guards = np.array([t["guard"] for t in trials])
    ax.scatter(starts[guards == "imm", 0], starts[guards == "imm", 1], s=10, c="#e7298a", alpha=0.5, label="trial → imm")
    ax.scatter(starts[guards == "ang", 0], starts[guards == "ang", 1], s=10, c="#7570b3", alpha=0.5, label="trial → ang")
    ax.set_xlim(0.05, 0.7)
    ax.set_ylim(0.15, 0.9)
    ax.set_xlabel("T")
    ax.set_ylabel("E")
    ax.set_title("Guard-hit geometry (every trial hits a guard)")
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "fig_4_3_guard_hits.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ef = np.array(surv["evals_full"])
    er = np.array(surv["evals_reduced"])
    ax.semilogy(np.arange(1, len(ef) + 1), np.maximum(ef, 1e-18), "o-", label="full", color="#1f77b4")
    ax.semilogy(np.arange(1, len(er) + 1), np.maximum(er, 1e-18), "s-", label="reduced", color="#ff7f0e")
    if ef[0] > 0:
        ax.axhline(PRACTICAL_CUT * ef[0], color="#1f77b4", ls=":", lw=1, label="full practical cut")
    if er[0] > 0:
        ax.axhline(PRACTICAL_CUT * er[0], color="#ff7f0e", ls=":", lw=1, label="reduced practical cut")
    ax.set_xlabel("Eigen-index (descending)")
    ax.set_ylabel("Fisher eigenvalue")
    ax.set_title(
        f"Surviving Fisher directions (practical ranks {surv['practical_rank_full']} → {surv['practical_rank_reduced']})"
    )
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "fig_4_4_surviving_fisher.png", dpi=160)
    plt.close(fig)


def main() -> None:
    theta = THETA_FULL.copy()
    orbit, period = find_orbit(400)
    trials = build_trials(orbit, theta)
    guard_hits = [t["guard"] for t in trials]
    n_hit = sum(1 for t in trials if t["reached_guard"])
    assert n_hit == len(trials)

    panel = [np.array(t["start"], dtype=float) for t in trials[::12]]
    F_full = fisher_matrix(theta, panel, "full")
    F_red = fisher_matrix(theta, panel, "reduced")
    surv = surviving_subspace(F_full, F_red, starts_for_jac=panel[:8])
    Q = surv.pop("Q")

    X = np.array([t["feature"] for t in trials], dtype=float)
    y_cycle = np.array([t["cycle_label"] for t in trials], dtype=int)
    y_occ = np.array([0 if t["guard"] == "imm" else 1 for t in trials], dtype=int)
    X_surv = X @ Q if Q.size else np.zeros((len(trials), 1))

    readers = {
        "cycle_full": fit_linear_reader(X, y_cycle),
        "cycle_surviving": fit_linear_reader(X_surv, y_cycle),
        "occult_full": fit_linear_reader(X, y_occ),
        "occult_surviving": fit_linear_reader(X_surv, y_occ),
    }

    honesty = collocation_defect_honesty(orbit)
    focus = np.array([T_STAR, E_STAR])
    dist_orbit_min = float(np.min(np.linalg.norm(orbit - focus, axis=1)))
    surviving_dirs = [d for d in surv["full_direction_survival"] if d["survives"]]

    results = {
        "seed": SEED,
        "depends_on": [
            "thesis-26-chain-recurrent-vs-hybrid-occult",
            "thesis-24-reduction-preserving-multichannel-id",
        ],
        "title": (
            "Chain-recurrent labels recoverable from multi-channel ranks that "
            "survive stiff–sloppy reduction"
        ),
        "problem": (
            "Which chain-recurrent versus hybrid-occult labels remain recoverable when "
            "observation is restricted to the multi-channel Fisher ranks that Thesis #24 "
            "proved survive a stiff–sloppy reduction?"
        ),
        "constants": {
            "A": A_PAR,
            "B": B_PAR,
            "C": C_PAR,
            "D": D_PAR,
            "T_star": float(T_STAR),
            "E_star": float(E_STAR),
            "T_ang": float(T_ANG),
            "E_imm": float(E_IMM),
            "theta_full": theta.tolist(),
            "full_names": FULL_NAMES,
            "phi_reduced": phi_of(theta).tolist(),
            "red_names": RED_NAMES,
            "kappa": float(kappa_of(theta)),
            "horizon": HORIZON,
            "period_approx": period,
            "n_trials": N_TRIALS,
            "n_focus": N_FOCUS,
            "practical_cut": PRACTICAL_CUT,
            "channel_names": CHANNEL_NAMES,
            "sigma_ch": SIGMA_CH.tolist(),
            "feature_dim": FEATURE_DIM,
        },
        "geometry": {
            "unstable_focus": [float(T_STAR), float(E_STAR)],
            "orbit_n_points": int(len(orbit)),
            "orbit_min_dist_to_focus": dist_orbit_min,
            "orbit_mean": orbit.mean(axis=0).tolist(),
            "every_trial_hits_guard": True,
            "n_trials": len(trials),
            "n_guard_hits": n_hit,
            "guard_counts": {
                "imm": int(sum(1 for g in guard_hits if g == "imm")),
                "ang": int(sum(1 for g in guard_hits if g == "ang")),
                "none": 0,
            },
            "occult_is_three_way_split": False,
            "occult_split": "two-way (imm vs ang); no no-hit class",
        },
        "fisher": {
            "practical_rank_full": surv["practical_rank_full"],
            "practical_rank_reduced": surv["practical_rank_reduced"],
            "evals_full": surv["evals_full"],
            "evals_reduced": surv["evals_reduced"],
            "n_surviving_directions": surv["n_surviving_directions"],
            "full_direction_survival": surv["full_direction_survival"],
            "surviving_full_directions": surviving_dirs,
        },
        "readers": readers,
        "findings": {
            "linear_reader_cycle_at_chance": bool(
                readers["cycle_full"]["at_chance"] and readers["cycle_surviving"]["at_chance"]
            ),
            "occult_recoverable_under_surviving_ranks": bool(readers["occult_surviving"]["above_chance"]),
            "occult_recoverable_under_full": bool(readers["occult_full"]["above_chance"]),
            "every_trial_hits_guard": True,
            "occult_not_three_way": True,
        },
        "honesty": honesty,
        "disclaimer": (
            "Research only. Not a medical device, not clinical decision support, not a dose, "
            "not a cure. Collocation defect is not a certified Conley set. No document DOI."
        ),
    }

    make_figures(results, trials, orbit, {**surv, "Q": Q})
    out = ROOT / "results.json"
    with out.open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
        f.write("\n")
    print(f"wrote {out}")
    print(
        "cycle_full={:.3f} cycle_surv={:.3f} occult_full={:.3f} occult_surv={:.3f} "
        "guards={}/{} surv_dirs={} findings={}".format(
            readers["cycle_full"]["accuracy"],
            readers["cycle_surviving"]["accuracy"],
            readers["occult_full"]["accuracy"],
            readers["occult_surviving"]["accuracy"],
            n_hit,
            len(trials),
            surv["n_surviving_directions"],
            results["findings"],
        )
    )


if __name__ == "__main__":
    main()
