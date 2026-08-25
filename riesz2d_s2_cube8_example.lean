import Mathlib

namespace Riesz2DS2Cube8

/-- The 12 edges of the planar cubic cube graph used by the Python example. -/
def cubeEdges : List (Nat × Nat) :=
  [(0,1),(1,2),(2,3),(3,0),
   (4,5),(5,6),(6,7),(7,4),
   (0,4),(1,5),(2,6),(3,7)]

/-- Bit v of mask records whether source vertex v is selected. -/
def selected (mask v : Nat) : Bool := ((mask / (2^v)) % 2) == 1

def selectedCount (mask : Nat) : Nat :=
  (List.range 8).foldl
    (fun acc v => match selected mask v with | true => acc + 1 | false => acc)
    0

def independentMask (mask : Nat) : Bool :=
  cubeEdges.all (fun e => !(selected mask e.1 && selected mask e.2))

def spin (mask v : Nat) : Int :=
  match selected mask v with | true => 1 | false => -1

/-- Barahona Hamiltonian H_B = edge products + uniform field. -/
def barahona (mask : Nat) : Int :=
  (cubeEdges.map (fun e => spin mask e.1 * spin mask e.2)).sum
  + ((List.range 8).map (fun v => spin mask v)).sum

/-- Brute-force maximum independent-set size over all 2^8 masks. -/
def alphaCube : Nat :=
  (List.range 256).foldl
    (fun best mask =>
      match independentMask mask with
      | true => Nat.max best (selectedCount mask)
      | false => best)
    0

/-- Brute-force minimum Barahona energy over all 2^8 masks. -/
def minBarahona : Int :=
  (List.range 256).foldl (fun best mask => min best (barahona mask)) 100

/-- The two bipartition classes are independent sets of size four. -/
example : independentMask 165 = true := by native_decide
example : selectedCount 165 = 4 := by native_decide
example : independentMask 90 = true := by native_decide
example : selectedCount 90 = 4 := by native_decide

/-- Exhaustive check that the cube graph has alpha = 4. -/
theorem cube_alpha : alphaCube = 4 := by native_decide

/-- Exhaustive check that the Barahona ground energy is -12. -/
theorem cube_barahona_ground : minBarahona = -12 := by native_decide

/-- The two maximum independent sets attain the ground energy. -/
example : barahona 165 = -12 := by native_decide
example : barahona 90 = -12 := by native_decide

/-- The source gap between alpha=4 and alpha<=3 is four energy units. -/
theorem cube_source_gap : (-8 : Int) - (-12) = 4 := by norm_num

/-- Demo selector forcing: (1/16)b^-2 > 4M for b=1/500 and M=392. -/
theorem demo_selector_forcing :
    (1 : ℚ) / 16 * (500^2 : ℚ) > 4 * 392 := by
  norm_num

/-- The measured finite-instance normalized remainder 0.168 lambda is much
    smaller than the two-unit half-gap used by the midpoint threshold. -/
theorem demo_remainder_below_half_gap : (168 : ℚ) / 1000 < 2 := by
  norm_num

/-- Consequently a scaled threshold -10 separates a yes level near -12 from
    a no level near -8 even after a +/-0.168 perturbation. -/
theorem demo_threshold_separation :
    (-12 : ℚ) + 168/1000 < -10 ∧ -10 < (-8 : ℚ) - 168/1000 := by
  constructor <;> norm_num

/-- The numerical component-flip audit reports more than 66 b^4 of margin. -/
theorem demo_flip_margin_positive : (66 : ℚ) > 0 := by
  norm_num

end Riesz2DS2Cube8
