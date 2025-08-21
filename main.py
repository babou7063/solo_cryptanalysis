import time
import matplotlib.pyplot as plt
import numpy as np
from elliptic_curve import Point, EllipticCurve
from basic_rho import pollard_rho
from additive_walk_rho import retry_walks, is_distinguished
from rho_negation_mapV2 import NegationMapRho

def find_curve_point(curve, start_x=1):
    """Find a valid point on the elliptic curve."""
    for x in range(start_x, curve.p):
        y_squared = (x**3 + curve.a * x + curve.b) % curve.p
        for y in range(curve.p):
            if (y * y) % curve.p == y_squared:
                return Point(x, y, curve)
    raise ValueError("No valid point found")

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
        'basic_rho_times': [],
        'additive_walk_times': [],
        'negation_map_times': [],
        'basic_rho_success': [],
        'additive_walk_success': [],
        'negation_map_success': []
    }
    
    print("=== Pollard's Rho Time Complexity Analysis ===\n")
    
    for p, a, b in test_cases:
        #print(f"Testing with prime p = {p}")
        curve = EllipticCurve(a, b, p)
        
        # Find a point on the curve
        P = find_curve_point(curve)
        order = curve.find_order(P)
        
        # Choose a secret k (not too large to ensure solvability)
        k_secret = min(7, order // 2)
        Q = curve.scalar_mul(k_secret, P)
        
        #print(f"  Point P: {P}")
        #print(f"  Order: {order}")
        #print(f"  Secret k: {k_secret}")
        #print(f"  Target Q: {Q}")
        
        results['primes'].append(p)
        
        # Test 1: Basic Rho
        print(" Testing Basic Rho...")
        start_time = time.time()
        try:
            k_found, _, _, steps = pollard_rho(P, Q, order, curve)
            end_time = time.time()
            
            basic_time = end_time - start_time
            basic_success = (k_found == k_secret)
            
            print(f"    Time: {basic_time:.4f}s, Steps: {steps}, Success: {basic_success}")
            results['basic_rho_times'].append(basic_time)
            results['basic_rho_success'].append(basic_success)
            
        except Exception as e:
            end_time = time.time()
            basic_time = end_time - start_time
            #print(f"    Failed: {e} (Time: {basic_time:.4f}s)")
            results['basic_rho_times'].append(basic_time)
            results['basic_rho_success'].append(False)
        
        # Test 2: Additive Walk Rho
        print("  Testing Additive Walk Rho...")
        start_time = time.time()
        try:
            k_found = retry_walks(P, Q, order, curve, r=8, is_distinguished=is_distinguished, max_attempts=5)
            end_time = time.time()
            
            additive_time = end_time - start_time
            additive_success = (k_found == k_secret) if k_found is not None else False
            
            print(f"    Time: {additive_time:.4f}s, Success: {additive_success}")
            results['additive_walk_times'].append(additive_time)
            results['additive_walk_success'].append(additive_success)
            
        except Exception as e:
            end_time = time.time()
            additive_time = end_time - start_time
            print(f"    Failed: {e} (Time: {additive_time:.4f}s)")
            results['additive_walk_times'].append(additive_time)
            results['additive_walk_success'].append(False)
        
        # Test 3: Negation Map Rho
        print(" Testing Negation Map Rho...")
        start_time = time.time()
        try:
            solver = NegationMapRho(curve, P, Q, order, r=32)
            k_found = solver.solve_ecdlp(max_walks=100)
            end_time = time.time()
            
            negation_time = end_time - start_time
            negation_success = (k_found == k_secret) if k_found is not None else False
            
            print(f"    Time: {negation_time:.4f}s, Success: {negation_success}")
            results['negation_map_times'].append(negation_time)
            results['negation_map_success'].append(negation_success)
            
        except Exception as e:
            end_time = time.time()
            negation_time = end_time - start_time
            print(f"    Failed: {e} (Time: {negation_time:.4f}s)")
            results['negation_map_times'].append(negation_time)
            results['negation_map_success'].append(False)
        
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

if __name__ == "__main__":
    
    results = run_time_complexity_analysis()
    
    #print_summary(results)
    
    plot_results(results)