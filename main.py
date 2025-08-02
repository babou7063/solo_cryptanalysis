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
        print(f"Testing with prime p = {p}")
        curve = EllipticCurve(a, b, p)
        
        # Find a point on the curve
        P = find_curve_point(curve)
        order = curve.find_order(P)
        
        # Choose a secret k (not too large to ensure solvability)
        k_secret = min(7, order // 2)
        Q = curve.scalar_mul(k_secret, P)
        
        print(f"  Point P: {P}")
        print(f"  Order: {order}")
        print(f"  Secret k: {k_secret}")
        print(f"  Target Q: {Q}")
        
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
            print(f"    Failed: {e} (Time: {basic_time:.4f}s)")
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
    
    print(f"  Basic Rho: {basic_success_rate:.1f}%")
    print(f"  Additive Walk: {additive_success_rate:.1f}%")
    print(f"  Negation Map: {negation_success_rate:.1f}%")
    
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
    
    print("\nAverage execution time (all runs):")
    if basic_avg_all: print(f"  Basic Rho: {basic_avg_all:.4f}s")
    if additive_avg_all: print(f"  Additive Walk: {additive_avg_all:.4f}s") 
    if negation_avg_all: print(f"  Negation Map: {negation_avg_all:.4f}s")
    
    print("\nAverage execution time (successful runs only):")
    if basic_avg: print(f"  Basic Rho: {basic_avg:.4f}s")
    if additive_avg: print(f"  Additive Walk: {additive_avg:.4f}s") 
    if negation_avg: print(f"  Negation Map: {negation_avg:.4f}s")


if __name__ == "__main__":
    
    results = run_time_complexity_analysis()
    
    print_summary(results)