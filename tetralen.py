import numpy as np
from numpy import sin, cos, sqrt, arcsin, arccos, pi, isfinite, array, dot, clip, sort
import time
import random


# Created by Mohammed Abdellateef
# Root-finding driver rewritten for correctness; the geometric primitives
# (midUp, midDown, thirdLen, riTri, isoTri) are unchanged/validated against
# known geometry. The coarse bracket scan is vectorized with NumPy for
# speed; final roots are refined with scalar bisection to `precision`.

def tetraLen(x1, x2, x3, ph1, ph2, ph3, precision, samples=8000):
   """
   Tetrahedreon lengths:
   	finding missing lengths in tetrahedreon from known tetrahedreon
      base lengths (3 lengths) and head-base angles (3 angles).

   In: 3 lengths, 3 angles
   Out: 3 lengths

   Convention: base triangle vertices Left (L), Mid (M), Right (R), apex H.
      x1 = |LM|, x3 = |MR|, x2 = |LR|
      ph1 = angle L-H-M, ph3 = angle M-H-R, ph2 = angle L-H-R
   Returns [HM, HL, HR] for each valid solution found.
   """
   T = time.time()

   # `precision` is accepted in two conventions for backward compatibility:
   #   - a decimal-digit count (e.g. 4, matching the original code's
   #     round(x, 4) usage) -- values >= 1 are treated this way
   #   - a raw absolute tolerance (e.g. 1e-4) -- values < 1 are used as-is
   # Everything below uses the normalized epsilon `tol`, never the raw
   # `precision` argument, so a caller passing "4" doesn't silently turn
   # into a tolerance of 4 units (which previously caused early bisection
   # termination and, far worse, made the solution-dedup radius huge
   # enough to merge genuinely different valid solutions into one).
   tol = 10 ** (-precision) if precision >= 1 else precision
   dedup_radius = max(1e-6, tol)

   # Farthest Points (longest lengths):
   # Equally sides Tetrahedreon
   def isoTri(x1, ph1):
      th1 = (pi-ph1)/2
      lr = x1*sin(th1)/sin(ph1)
      return lr

   # longest length Projection
   def riTri(x1, ph1):
      inside = x1/sin(ph1)     # hypotenuse -> upper bound on the shared side (HM)
      thp1 = pi/2-ph1
      outside = inside*sin(thp1)
      return outside, inside

   # Scalar versions (used for the final bisection refinement).
   def midUp(l, x, ph):
      ath = l*sin(ph)/x
      if ath > 1:
         return np.nan
      th = arcsin(ath)
      if pi-th > th:
         th = pi - th
      thu = th - ph
      return l*sin(thu)/sin(th)

   def midDown(l, x, ph):
      ath = l*sin(ph)/x
      if ath > 1:
         return np.nan
      th = arcsin(ath)
      if pi-th > th:
         th = pi - th
      thl = th + ph
      return l*sin(thl)/sin(th)

   # Vectorized versions (used for the fast coarse scan).
   def midUp_vec(l, x, ph):
      with np.errstate(divide='ignore', invalid='ignore'):
        ath = l * sin(ph) / x
        valid = np.isfinite(ath) & (np.abs(ath) <= 1)

        ath_c = np.clip(ath, -1.0, 1.0)
        th = arcsin(ath_c)
        th = np.where(pi - th > th, pi - th, th)

        thu = th - ph
        den = sin(th)

        out = np.divide(
            l * sin(thu),
            den,
            out=np.full_like(den, np.nan, dtype=float),
            where=np.abs(den) > 1e-14
        )

        return np.where(valid, out, np.nan)


   def midDown_vec(l, x, ph):
      with np.errstate(divide='ignore', invalid='ignore'):
        ath = l * sin(ph) / x
        valid = np.isfinite(ath) & (np.abs(ath) <= 1)

        ath_c = np.clip(ath, -1.0, 1.0)
        th = arcsin(ath_c)
        th = np.where(pi - th > th, pi - th, th)

        thl = th + ph
        den = sin(th)

        out = np.divide(
            l * sin(thl),
            den,
            out=np.full_like(den, np.nan, dtype=float),
            where=np.abs(den) > 1e-14
        )

        return np.where(valid, out, np.nan)
        
   # cos rule to find missing length
   def thirdLen(l1, l2, th):
      return sqrt(l1**2 + l2**2 - 2*l1*l2*cos(th))

   # --- upper bound on HM (the shared/searched side) ---
   # Each of the two adjacent triangles (L-H-M via x1,ph1 and M-H-R via
   # x3,ph3) independently caps how large HM can be; the true cap is
   # whichever is smaller.
   _, inl = riTri(x1, ph1)   # inl = x1/sin(ph1): HM cap from the L-H-M triangle
   _, inr = riTri(x3, ph3)   # inr = x3/sin(ph3): HM cap from the M-H-R triangle
   hm_max = min(inl, inr)

   sol = []
   if not (isfinite(hm_max) and hm_max > 0):
      elapsed = time.time() - T
      return sol, 0, elapsed

   # --- scan the 4 branches (HL from midUp/midDown x HR from midUp/midDown) ---
   # For each branch, g(HM) = thirdLen(HL(HM), HR(HM), ph2) - x2 is continuous
   # over (0, hm_max]; bracket sign changes (vectorized) then bisect each
   # bracket (scalar) to `precision`.
   branches = [(midUp, midUp, midUp_vec, midUp_vec),
               (midUp, midDown, midUp_vec, midDown_vec),
               (midDown, midUp, midDown_vec, midUp_vec),
               (midDown, midDown, midDown_vec, midDown_vec)]
   eps = hm_max * 1e-9
   hs = np.linspace(eps, hm_max, samples)

   def g_scalar(fL, fR, h):
      hl = fL(h, x1, ph1)
      hr = fR(h, x3, ph3)
      if not (isfinite(hl) and isfinite(hr)):
         return np.nan, None, None
      return thirdLen(hl, hr, ph2) - x2, hl, hr

   found = []
   for fL, fR, fL_vec, fR_vec in branches:
      hl_arr = fL_vec(hs, x1, ph1)
      hr_arr = fR_vec(hs, x3, ph3)
      valid = isfinite(hl_arr) & isfinite(hr_arr)
      v_arr = np.where(valid, thirdLen(hl_arr, hr_arr, ph2) - x2, np.nan)

      v0s, v1s = v_arr[:-1], v_arr[1:]
      both_finite = isfinite(v0s) & isfinite(v1s)
      exact = both_finite & (v0s == 0)
      crossing = both_finite & (v0s * v1s < 0)
      bracket_idx = np.nonzero(exact | crossing)[0]

      for i in bracket_idx:
         if exact[i]:
            root = hs[i]
         else:
            a, b, fa = hs[i], hs[i+1], v0s[i]
            for _ in range(200):
               m = (a + b) / 2.0
               fm, _, _ = g_scalar(fL, fR, m)
               if not isfinite(fm):
                  break
               if (fa < 0) == (fm < 0):
                  a, fa = m, fm
               else:
                  b = m
               if (b - a) < 1e-13 * max(1.0, hm_max):
                  break
            root = (a + b) / 2.0
         _, hl_r, hr_r = g_scalar(fL, fR, root)
         if hl_r is not None and hr_r is not None and min(root, hl_r, hr_r) > 0:
            found.append([root, hl_r, hr_r])

   # dedupe near-identical solutions (different branches can converge to
   # the same physical point, e.g. at hm_max where up == down)
   for cand in found:
      if not any(all(abs(a - b) < dedup_radius for a, b in zip(cand, existing))
                 for existing in sol):
         sol.append(cand)

   elapsed = time.time() - T
   nSol = len(sol)
   return sol, nSol, elapsed
   
# ============================================================
# TETRALEN TEST DATA GENERATOR
# ============================================================

def generate_tetra_test():
    """
    Generate a random tetrahedron and convert it to the
    input format required by tetraLen().

    Returns:
        x1, x2, x3          known base lengths
        ph1, ph2, ph3       head-base angles
        expected            three missing lengths
        points              original 3-D coordinates
    """

    # Random tetrahedron vertices
    A = array([0.0, 0.0, 0.0])

    B = array([
        random.uniform(2.0, 10.0),
        0.0,
        0.0
    ])

    C = array([
        random.uniform(-5.0, 10.0),
        random.uniform(2.0, 10.0),
        0.0
    ])

    D = array([
        random.uniform(-5.0, 10.0),
        random.uniform(-5.0, 10.0),
        random.uniform(2.0, 10.0)
    ])

    # Distance helper
    def dist(P, Q):
        return sqrt(sum((P - Q) ** 2))

    # --------------------------------------------------------
    # Six actual tetrahedron edges
    #
    # Base:
    #   AB = x1
    #   AC = x2
    #   BC = x3
    #
    # Unknown:
    #   AD
    #   BD
    #   CD
    # --------------------------------------------------------

    AB = dist(A, B)
    AC = dist(A, C)
    BC = dist(B, C)

    AD = dist(A, D)
    BD = dist(B, D)
    CD = dist(C, D)

    # --------------------------------------------------------
    # Angles at the head D
    #
    # ph1 = angle A-D-B
    # ph2 = angle A-D-C
    # ph3 = angle B-D-C
    # --------------------------------------------------------

    def angle(P, Q, R):
        """
        Angle P-Q-R.
        """
        v1 = P - Q
        v2 = R - Q

        c = dot(v1, v2) / (dist(P, Q) * dist(R, Q))
        

        # protect against floating-point errors
        c = clip(c, -1.0, 1.0)

        return arccos(c)

    ph1 = angle(A, D, B)
    ph2 = angle(A, D, C)
    ph3 = angle(B, D, C)

    return (
        AB, AC, BC,
        ph1, ph2, ph3,
        array([AD, BD, CD]),
        (A, B, C, D)
    )


def test_tetraLen(n=10, precision=4):
    """
    Generate and test n random tetrahedra against tetraLen().
    """

    print("=" * 70)
    print("TETRALEN RANDOM TEST")
    print("=" * 70)

    passed = 0
    failed = 0

    for test_no in range(1, n + 1):

        (
            x1, x2, x3,
            ph1, ph2, ph3,
            expected,
            points
        ) = generate_tetra_test()

        print("\nTest", test_no)
        print("-" * 70)

        print("Known base:")
        print("x1 =", x1)
        print("x2 =", x2)
        print("x3 =", x3)

        print("\nAngles:")
        print("ph1 =", ph1)
        print("ph2 =", ph2)
        print("ph3 =", ph3)

        print("\nExpected missing lengths:")
        print("AD =", expected[0])
        print("BD =", expected[1])
        print("CD =", expected[2])

        try:
            result, nsol, elapsed = tetraLen(
                x1, x2, x3,
                ph1, ph2, ph3,
                precision
            )

            print("\nTetralen result:")
            print(result)
            print("Number of solutions:", nsol)
            print("Time:", elapsed)

            # Check whether any returned solution matches
            # the actual tetrahedron.
            found = False

            for sol in result:

                if len(sol) != 3:
                    continue

                # Because tetraLen may return a different
                # ordering depending on its geometric branch,
                # compare sorted lengths.
                a = sort(array(sol, dtype=float))
                b = sort(expected)

                if all(abs(a - b) < 10 ** (-precision)):
                    found = True
                    break

            if found:
                print("PASS")
                passed += 1
            else:
                print("FAIL")
                print("Expected:", expected)
                print("Returned:", result)
                failed += 1

        except Exception as e:
            print("\nERROR:", type(e).__name__, e)
            failed += 1

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("Passed:", passed)
    print("Failed:", failed)
    print("Total :", n)

if __name__ == "__main__":
    test_tetraLen(10)
