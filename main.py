import time
import matplotlib.pyplot as plt
import numpy as np
import random
import math
from elliptic_curve import Point, EllipticCurve
from basic_rho import pollard_rho
from additive_walk_rho import retry_walks, is_distinguished
from rho_negation_mapV2 import NegationMapRho

def find_point_with_min_order(curve, min_order=None, max_try=5_000):
    """
    Find a point on the elliptic curve with an order greater or equal to a given minimum.
    """
    if min_order is None:
        min_order = max(50, int(curve.p ** 0.5))

    tries = 0
    for x in range(curve.p):
        y2 = (x**3 + curve.a * x + curve.b) % curve.p
        for y in range(curve.p):
            if (y * y) % curve.p == y2:
                P = Point(x, y, curve)
                try:
                    order = curve.find_order(P, max_iter=10_000)
                except Exception:
                    continue
                if order >= min_order:
                    return P, order
                tries += 1
                if tries >= max_try:
                    raise ValueError("Impossible to find a point with sufficient order.")
    raise ValueError("No valid point found on the curve.")

def largest_prime_factor(n: int) -> int:
    f, lp = 2, 1
    while f * f <= n:
        while n % f == 0:
            lp, n = f, n // f
        f = 3 if f == 2 else f + 2
    return n if n > 1 else lp

def robust_basic_rho(P, Q, order, curve, r=3, retries=8):
    """
    A robust version of the basic Pollard's rho algorithm that
    runs the algorithm multiple times and returns the first successful
    result. If all attempts fail, return None, None, False, and a dictionary
    containing the number of operations performed.
    """
    last_err = None
    for _ in range(retries):
        try:
            k_found, _, _, steps, ops = pollard_rho(P, Q, order, curve, r=r)
            return k_found, steps, True, ops
        except Exception as e:
            last_err = e
            continue
    return None, None, False, {'adds':0,'dbls':0,'total':0}

def make_is_distinguished(n, alpha=8, min_bits=1):
    t_bits = max(min_bits, int(round(math.log2(max(2, int(math.isqrt(n) // max(1, alpha)))))))
    mask = (1 << t_bits) - 1
    def _pred(W):
        return (not W.at_infinity) and ((W.x & mask) == 0)
    return _pred

def run_time_complexity_analysis():
    """Run time complexity analysis for all three algorithms."""
    
    # Test parameters - start small and increase
    test_cases = [
        (101, 2, 3),    # Small prime
        (211, 2, 3),    # Medium prime  
        (307, 2, 3),    # Larger prime
        (401, 2, 3),    # Even larger
        (503, 2, 3),    # Larger still
    ]
    
    results = {
        'primes': [],
        'orders': [],
        'basic_rho_times': [],
        'additive_walk_times': [],
        'negation_map_times': [],
        'basic_rho_success': [],
        'additive_walk_success': [],
        'negation_map_success': [],
        'basic_ops': [], 
        'additive_ops': [], 
        'negation_ops': [],
        'basic_ops_precompute': [], 
        'additive_ops_precompute': [], 
        'negation_ops_precompute': []
    }
    
    print("=== Pollard's Rho Time Complexity Analysis ===\n")
    
    for p, a, b in test_cases:
        curve = EllipticCurve(a, b, p)

        P, order = find_point_with_min_order(curve)
        ell = largest_prime_factor(order)
        if ell >= max(int(math.sqrt(order)), 16):
            h = order // ell
            P = curve.scalar_mul(h, P)
            order = ell
        k_secret = random.randrange(1, order)
        Q = curve.scalar_mul(k_secret, P)

        results['primes'].append(p)
        results['orders'].append(order)
        print(f"[p={p}] ordre(P)={order}, k_secret={k_secret}")

        # --- Basic Rho (robuste)
        print("  Testing Basic Rho (robust)...")
        t0 = time.time()
        k_found, steps, ok, ops = robust_basic_rho(P, Q, order, curve, r=3, retries=10)
        t1 = time.time()
        basic_time = t1 - t0
        basic_success = (ok and ((k_found - k_secret) % order == 0))
        print(f"    Time: {basic_time:.4f}s, Steps: {steps}, Success: {basic_success}")
        results['basic_rho_times'].append(basic_time)
        results['basic_rho_success'].append(basic_success)
        if basic_success:
            results['basic_ops'].append(ops['total'])
            results['basic_ops_precompute'].append(0)  # pas de précomp pour basic
        else:
            results['basic_ops'].append(None); results['basic_ops_precompute'].append(None)


        # --- Additive Walk
        print("  Testing Additive Walk Rho...")
        t0 = time.time()
        try:
            dp_pred = make_is_distinguished(order, alpha=2, min_bits=2)
            k_found, stats = retry_walks(P, Q, order, curve, r=8, is_distinguished=dp_pred, max_attempts=6, return_stats=True)
            additive_success = (k_found is not None and ((k_found - k_secret) % order == 0))
            if additive_success:
                results['additive_ops'].append(stats['ops_total'])
                results['additive_ops_precompute'].append(stats['ops_precompute'])
            else:
                results['additive_ops'].append(None); results['additive_ops_precompute'].append(None)
        except Exception:
            additive_success = False
        t1 = time.time()
        additive_time = t1 - t0
        print(f"    Time: {additive_time:.4f}s, Success: {additive_success}")
        results['additive_walk_times'].append(additive_time)
        results['additive_walk_success'].append(additive_success)

        # --- Negation Map
        print("  Testing Negation Map Rho...")
        t0 = time.time()
        try:
            solver = NegationMapRho(curve, P, Q, order, r=32)
            k_found, stats = solver.solve_ecdlp(max_walks=800, return_stats=True)
            negation_success = (k_found is not None and ((k_found - k_secret) % order == 0))
            if negation_success:
                results['negation_ops'].append(stats['ops_total'])
                results['negation_ops_precompute'].append(stats['ops_precompute'])
            else:
                results['negation_ops'].append(None); results['negation_ops_precompute'].append(None)

        except Exception:
            negation_success = False
        t1 = time.time()
        negation_time = t1 - t0
        print(f"    Time: {negation_time:.4f}s, Success: {negation_success}")
        results['negation_map_times'].append(negation_time)
        results['negation_map_success'].append(negation_success)

        print()

    return results

def print_summary(results):
    """Print a summary of the results."""
    print("=== SUMMARY ===")
    print(f"Test cases: {len(results['primes'])}")
    print(f"Primes tested: {results['primes']}")
    
    print("\nSuccess rates:")
    basic_success_rate = sum(results['basic_rho_success']) / len(results['basic_rho_success']) * 100
    additive_success_rate = sum(results['additive_walk_success']) / len(results['additive_walk_success']) * 100
    negation_success_rate = sum(results['negation_map_success']) / len(results['negation_map_success']) * 100
    
    print(f"==> Basic Rho: {basic_success_rate:.1f}%")
    print(f"==> Additive Walk: {additive_success_rate:.1f}%")
    print(f"==> Negation Map: {negation_success_rate:.1f}%")
    
    # Calculate average times for all runs (successful and failed)
    def avg_time(times):
        valid_times = [t for t in times if t is not None]
        return np.mean(valid_times) if valid_times else None
    
    # Calculate average times for successful runs only
    def avg_successful_time(times, successes):
        successful_times = [t for t, s in zip(times, successes) if t is not None and s]
        return np.mean(successful_times) if successful_times else None
    
    basic_avg_all = avg_time(results['basic_rho_times'])
    additive_avg_all = avg_time(results['additive_walk_times'])
    negation_avg_all = avg_time(results['negation_map_times'])
    
    basic_avg = avg_successful_time(results['basic_rho_times'], results['basic_rho_success'])
    additive_avg = avg_successful_time(results['additive_walk_times'], results['additive_walk_success'])
    negation_avg = avg_successful_time(results['negation_map_times'], results['negation_map_success'])
    
    print("\nAverage execution time (successful runs only):")
    if basic_avg: print(f"  Basic Rho: {basic_avg:.4f}s")
    if additive_avg: print(f"  Additive Walk: {additive_avg:.4f}s") 
    if negation_avg: print(f"  Negation Map: {negation_avg:.4f}s")
    
    print("\nAverage execution time (all runs):")
    if basic_avg_all: print(f"  Basic Rho: {basic_avg_all:.4f}s")
    if additive_avg_all: print(f"  Additive Walk: {additive_avg_all:.4f}s") 
    if negation_avg_all: print(f"  Negation Map: {negation_avg_all:.4f}s")



def plot_results(results):
    """Create plots for the time complexity analysis."""
    
    # Filter out None values for plotting
    def filter_data(times, primes):
        filtered_times = []
        filtered_primes = []
        for i, t in enumerate(times):
            if t is not None:
                filtered_times.append(t)
                filtered_primes.append(primes[i])
        return filtered_times, filtered_primes
    
    basic_times, basic_primes = filter_data(results['basic_rho_times'], results['primes'])
    additive_times, additive_primes = filter_data(results['additive_walk_times'], results['primes'])
    negation_times, negation_primes = filter_data(results['negation_map_times'], results['primes'])
    
    # Create subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Execution times comparison
    ax1.plot(basic_primes, basic_times, 'bo-', label='Basic Rho', linewidth=2, markersize=8)
    ax1.plot(additive_primes, additive_times, 'rs-', label='Additive Walk Rho', linewidth=2, markersize=8)
    ax1.plot(negation_primes, negation_times, 'g^-', label='Negation Map Rho', linewidth=2, markersize=8)
    ax1.set_xlabel('Prime p')
    ax1.set_ylabel('Execution Time (seconds)')
    ax1.set_title('Execution Time vs Prime Size')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')  # Log scale for better visualization
    
    # Plot 2: Success rates
    success_data = {
        'Basic Rho': sum(results['basic_rho_success']) / len(results['basic_rho_success']) * 100,
        'Additive Walk': sum(results['additive_walk_success']) / len(results['additive_walk_success']) * 100,
        'Negation Map': sum(results['negation_map_success']) / len(results['negation_map_success']) * 100
    }
    
    bars = ax2.bar(success_data.keys(), success_data.values(), color=['blue', 'red', 'green'], alpha=0.7)
    ax2.set_ylabel('Success Rate (%)')
    ax2.set_title('Algorithm Success Rates')
    ax2.set_ylim(0, 100)
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{height:.1f}%', ha='center', va='bottom')
    
    # Plot 3: Theoretical vs Actual complexity
    theoretical_times = [np.sqrt(np.pi * p / 4) / 1000 for p in results['primes']]  # Scaled theoretical
    
    ax3.plot(results['primes'], theoretical_times, 'k--', label='Theoretical O(√p)', linewidth=2)
    if basic_times:
        ax3.plot(basic_primes, basic_times, 'bo-', label='Basic Rho (Actual)', linewidth=2, markersize=6)
    if additive_times:
        ax3.plot(additive_primes, additive_times, 'rs-', label='Additive Walk (Actual)', linewidth=2, markersize=6)
    if negation_times:
        ax3.plot(negation_primes, negation_times, 'g^-', label='Negation Map (Actual)', linewidth=2, markersize=6)
    
    ax3.set_xlabel('Prime p')
    ax3.set_ylabel('Time (seconds)')
    ax3.set_title('Theoretical vs Actual Complexity')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_yscale('log')
    
    # Plot 4: Average execution time comparison
    avg_times = []
    labels = []
    colors = []
    
    if basic_times:
        avg_times.append(np.mean(basic_times))
        labels.append('Basic Rho')
        colors.append('blue')
    
    if additive_times:
        avg_times.append(np.mean(additive_times))
        labels.append('Additive Walk')
        colors.append('red')
    
    if negation_times:
        avg_times.append(np.mean(negation_times))
        labels.append('Negation Map')
        colors.append('green')
    
    if avg_times:
        bars = ax4.bar(labels, avg_times, color=colors, alpha=0.7)
        ax4.set_ylabel('Average Time (seconds)')
        ax4.set_title('Average Execution Time Comparison')
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax4.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                    f'{height:.4f}s', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.show()

def plot_ops(results):
    orders = results.get('orders', [])  # n effectifs

    def filt(xn, total, pre):
        """Filtre (n, total_ops, precomp_ops) -> (x triés, total, walk, pre)."""
        xs, tot, walk, prec = [], [], [], []
        for n, t, p in zip(xn, total, pre):
            if n is None or t is None or p is None:
                continue
            xs.append(int(n))
            tot.append(int(t))
            prec.append(int(p))
            walk.append(max(int(t) - int(p), 0))
        # tri par abscisse
        pairs = sorted(zip(xs, tot, walk, prec), key=lambda z: z[0])
        if not pairs:
            return [], [], [], []
        xs, tot, walk, prec = map(list, zip(*pairs))
        return xs, tot, walk, prec

    # Récupération et séparation walk / précomp (déjà triés par n)
    b_x, b_tot, b_walk, b_pre = filt(
        orders, results.get('basic_ops', []), results.get('basic_ops_precompute', [])
    )
    a_x, a_tot, a_walk, a_pre = filt(
        orders, results.get('additive_ops', []), results.get('additive_ops_precompute', [])
    )
    n_x, n_tot, n_walk, n_pre = filt(
        orders, results.get('negation_ops', []), results.get('negation_ops_precompute', [])
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))

    # --- Gauche : ops WALK vs n, avec réf c·√n ---
    if b_walk: ax1.plot(b_x, b_walk, 'bo-', label='Basic (walk)')
    if a_walk: ax1.plot(a_x, a_walk, 'rs-', label='Additive (walk)')
    if n_walk: ax1.plot(n_x, n_walk, 'g^-', label='Negation (walk)')

    # Référence c·√n calibrée sur la médiane(ops/√n) tous algos confondus
    xs_ref = sorted(set(b_x + a_x + n_x)) or sorted(set(int(x) for x in orders))
    if xs_ref:
        ratios = []
        for xs, ws in ((b_x, b_walk), (a_x, a_walk), (n_x, n_walk)):
            ratios += [w / math.sqrt(x) for x, w in zip(xs, ws) if x > 0 and w > 0]
        c = float(np.median(ratios)) if ratios else 1.0
        ax1.plot(xs_ref, [c * math.sqrt(x) for x in xs_ref], 'k--', label='c·√n')

    ax1.set_xscale('log'); ax1.set_yscale('log')
    ax1.set_xlabel('Subgroup order n')
    ax1.set_ylabel('Group operations (walk only)')
    ax1.set_title('Group ops per success (WALK only) vs n')
    ax1.grid(True, alpha=.3); ax1.legend()

    # --- Droite : barres empilées (moyenne PRECOMP + WALK) par algo ---
    def avg(arr):
        arr = [x for x in arr if x is not None]
        return float(np.mean(arr)) if arr else 0.0

    bars_labels, pre_vals, walk_vals = [], [], []

    if b_tot:
        bars_labels.append('Basic');    pre_vals.append(avg(b_pre)); walk_vals.append(avg(b_walk))
    if a_tot:
        bars_labels.append('Additive'); pre_vals.append(avg(a_pre)); walk_vals.append(avg(a_walk))
    if n_tot:
        bars_labels.append('Negation'); pre_vals.append(avg(n_pre)); walk_vals.append(avg(n_walk))

    x_pos = np.arange(len(bars_labels))
    ax2.bar(x_pos, pre_vals, label='Precompute', alpha=0.7)
    ax2.bar(x_pos, walk_vals, bottom=pre_vals, label='Walk', alpha=0.9)
    ax2.set_xticks(x_pos, bars_labels)
    ax2.set_ylabel('Avg group operations')
    ax2.set_title('Average ops per success (stacked: precomp + walk)')
    ax2.grid(True, axis='y', alpha=.3)
    ax2.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    
    results = run_time_complexity_analysis()
    
    #print_summary(results)
    plot_results(results)
    plot_ops(results)
    