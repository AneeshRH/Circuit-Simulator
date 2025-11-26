# graphz.py - Symbolic Laplace Solver (Exact Analytical Solution)
import numpy as np
import sympy as sym
import sys

# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================


def format_value(value, comp_type):
    """Format component values for display (e.g., 1000 -> 1.0kΩ)"""
    try:
        val = float(value)
        if comp_type == 'R':
            return f"{val/1e3:.1f}kΩ" if val >= 1e3 else f"{val:.1f}Ω"
        elif comp_type == 'C':
            return f"{val*1e6:.1f}μF" if val < 1e-3 else f"{val:.3f}F"
        elif comp_type == 'L':
            return f"{val*1e3:.1f}mH" if val < 1 else f"{val:.3f}H"
        return f"{val:.1f}" + ("V" if comp_type == 'V' else "A" if comp_type == 'I' else "")
    except:
        return str(value)


def parse_netlist_structure(netlist_str):
    """Used by UI to list components"""
    components = []
    for line in netlist_str.strip().split('\n'):
        parts = line.split()
        if len(parts) >= 4:
            components.append({
                'name': parts[0],
                'val_str': format_value(parts[3], parts[0][0]),
                'type': parts[0][0]
            })
    return components


def parse_value_str(val_str):
    """Parses engineering notation (1k, 10u, etc.)"""
    suffixes = {'T': 1e12, 'G': 1e9, 'M': 1e6, 'k': 1e3, 'm': 1e-3,
                'u': 1e-6, 'n': 1e-9, 'p': 1e-12, 'f': 1e-15}
    s = val_str.lower().rstrip('vahoHzfmΩ')
    try:
        if s[-1] in suffixes:
            mult = 1.0
            for k, v in suffixes.items():
                if s.endswith(k.lower()):
                    mult = v
                    s = s[:-1]
                    break
            return float(s) * mult
        return float(s)
    except:
        return 0.0

# ==============================================================================
# 2. SYMBOLIC LAPLACE SOLVER
# ==============================================================================


def solve_circuit_mna(netlist_str, t_start=0, t_end=0.1):
    """
    Solves the circuit using Modified Nodal Analysis in the LAPLACE DOMAIN.
    1. Converts components to Impedance (Z).
    2. Solves linear equations in 's'.
    3. Performs Inverse Laplace Transform to get exact f(t).
    """
    try:
        # --- A. PARSE NETLIST ---
        components = []
        nodes = set()

        for line in netlist_str.strip().split('\n'):
            parts = line.split()
            if len(parts) < 4:
                continue

            name, n1, n2, val_str = parts[0], int(
                parts[1]), int(parts[2]), parts[3]
            val = parse_value_str(val_str)
            nodes.update([n1, n2])

            components.append({
                'name': name, 'type': name[0].upper(),
                'n1': n1, 'n2': n2, 'val': val
            })

        if not components:
            return "Error: Empty netlist"

        # MNA Matrix Size
        num_nodes = max(nodes)  # Node 0 is ground
        m_size = num_nodes  # V1...Vn

        # Voltage sources add extra rows/cols to matrix
        v_sources = [c for c in components if c['type'] == 'V']
        num_vs = len(v_sources)
        dim = m_size + num_vs

        # --- B. BUILD MNA MATRIX IN S-DOMAIN ---
        s, t = sym.symbols('s t')

        # G (Conductance/Admittance) Matrix and I (Source) Vector
        G_mat = sym.zeros(dim, dim)
        I_vec = sym.zeros(dim, 1)

        # 1. Fill Passive Elements (Z)
        for comp in components:
            # 0-indexed for matrix (Node 1 is idx 0)
            n1, n2 = comp['n1'] - 1, comp['n2'] - 1
            val = comp['val']

            # Determine Admittance Y(s) = 1/Z(s)
            Y = 0
            if comp['type'] == 'R':
                Y = 1 / val
            elif comp['type'] == 'L':
                Y = 1 / (s * val)
            elif comp['type'] == 'C':
                Y = s * val
            elif comp['type'] == 'I':
                # Current source enters RHS vector
                # Step input: I(s) = val / s
                i_s = val / s
                if n1 >= 0:
                    I_vec[n1] -= i_s
                if n2 >= 0:
                    I_vec[n2] += i_s
                continue
            elif comp['type'] == 'V':
                continue  # Handle V sources separately

            # Stamp Y into Matrix
            if n1 >= 0:
                G_mat[n1, n1] += Y
            if n2 >= 0:
                G_mat[n2, n2] += Y
            if n1 >= 0 and n2 >= 0:
                G_mat[n1, n2] -= Y
                G_mat[n2, n1] -= Y

        # 2. Fill Voltage Sources
        # V sources add a constraint row: V_n1 - V_n2 = V_src
        for idx, vs in enumerate(v_sources):
            row = m_size + idx  # Extra row index
            n1, n2 = vs['n1'] - 1, vs['n2'] - 1
            val = vs['val']

            # Step Input: V(s) = val / s
            v_s = val / s

            # Stamp Incidence (1, -1)
            if n1 >= 0:
                G_mat[n1, row] = 1
                G_mat[row, n1] = 1
            if n2 >= 0:
                G_mat[n2, row] = -1
                G_mat[row, n2] = -1

            I_vec[row] = v_s

        # --- C. SOLVE SYMBOLICALLY ---
        # System: [G] * [X] = [I]
        # X = [V1(s), ..., Vn(s), I_Vsrc1(s)...]

        try:
            solution_s = G_mat.LUsolve(I_vec)
        except sym.matrices.common.NonInvertibleMatrixError:
            return "Error: Singular Matrix (Check Ground/Floating Nodes)"

        # --- D. INVERSE LAPLACE & EVALUATE ---
        t_points = np.linspace(t_start, t_end, 200)
        results = []

        # Helper to convert sym expr to numpy array
        def get_time_data(expr_s):
            # 1. Inverse Laplace
            expr_t = sym.inverse_laplace_transform(expr_s, s, t)
            # 2. Simplify (removes Heaviside(t) if t>0 assumed)
            # We assume t >= 0
            expr_t = expr_t.replace(sym.Heaviside(t), 1)
            # 3. Lambdify
            func_t = sym.lambdify(t, expr_t, modules=['numpy'])

            # Calculate
            vals = func_t(t_points)

            # Handle scalar result (if constant)
            if np.isscalar(vals):
                vals = np.full_like(t_points, vals)
            return vals

        # 1. Extract Node Voltages
        node_voltages = {}  # Map node_idx (1-based) to time array
        node_voltages[0] = np.zeros_like(t_points)  # Ground

        for i in range(num_nodes):
            v_s = solution_s[i]
            v_t_data = get_time_data(v_s)
            node_voltages[i+1] = v_t_data

            results.append({
                'name': f'v_{i+1}', 'val_str': 'Node', 'time': t_points,
                'voltage': v_t_data, 'current': np.zeros_like(t_points)
            })

        # 2. Extract Source Currents (from solution vector)
        src_currents = {}
        for idx, vs in enumerate(v_sources):
            row = m_size + idx
            i_s = solution_s[row]
            i_t_data = get_time_data(i_s)
            src_currents[vs['name']] = i_t_data

        # 3. Calculate Component Currents/Voltages
        for comp in components:
            n1, n2 = comp['n1'], comp['n2']
            v_drop = node_voltages[n1] - node_voltages[n2]

            curr = np.zeros_like(t_points)

            if comp['type'] == 'R':
                curr = v_drop / comp['val']
            elif comp['type'] == 'L':
                # I(s) = V(s) / sL
                # We can calculate this using the symbolic V_drop(s)
                # But we already have V(t), we can just use inverse laplace on V_drop(s)/sL for precision
                # Or numerically integrate. Let's do symbolic for consistency.

                # Re-calculate symbolic V_drop(s)
                v1_s = solution_s[n1-1] if n1 > 0 else 0
                v2_s = solution_s[n2-1] if n2 > 0 else 0
                v_drop_s = v1_s - v2_s

                i_L_s = v_drop_s / (s * comp['val'])
                curr = get_time_data(i_L_s)

            elif comp['type'] == 'C':
                # I(s) = V(s) * sC
                v1_s = solution_s[n1-1] if n1 > 0 else 0
                v2_s = solution_s[n2-1] if n2 > 0 else 0
                v_drop_s = v1_s - v2_s

                i_C_s = v_drop_s * (s * comp['val'])
                curr = get_time_data(i_C_s)

            elif comp['type'] == 'V':
                curr = src_currents[comp['name']]
            elif comp['type'] == 'I':
                curr = np.full_like(t_points, comp['val'])

            results.append({
                'name': comp['name'],
                'val_str': format_value(comp['val'], comp['type']),
                'time': t_points,
                'voltage': v_drop,
                'current': curr
            })

        return results

    except Exception as e:
        # Fallback for very complex circuits where symbolic fails
        return f"Symbolic Solver Error: {str(e)}"
