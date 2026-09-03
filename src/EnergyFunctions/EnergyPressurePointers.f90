!===================================================================
! Module: EnergyPressurePointers
! Purpose:
!   Provides unified interfaces to select appropriate subroutines
!   for energy evaluation, pressure calculation, and Rosenbluth
!   sampling based on the selected force field in the simulation.
!
!   Supported Force Fields and Their Mathematical Forms:
!
!   1. LJ_Q: Lennard-Jones and Electrostatic Potential
!      - U_LJ(r)           = ε * [ (σ / r)**12 − (σ / r)**6 ]
!      - U_Electrostatic   = (q1 * q2) / r
!
!   2. Mpipi: Wang–Frankel + Debye-Hückel Electrostatics
!      - U_WF(r)           = ε * α * [ ( (σ / r)**(2μ) − 1 ) − ( (3σ / r)**(2μ) − 1 ) ]**(2ν)
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!
!   3. HPS_single (HPS-Nucl): Hydropathy Scale Lennard-Jones + Debye-Hückel
!      - U_HPS(r)          = 4 * ε * λ * [ (σ / r)**12 − 2 * (σ / r)**6 ]
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!
!   4. HPS_piecewise (HPS-KR, KH, FB, TSCL-M2, Urry: Piecewise HPS LJ + Debye-Hückel
!      - u_LJ(r)           = 4 * ε * [ (σ / r)**12 − (σ / r)**6 ]
!      - U_HPS(r) =
!            u_LJ(r) + ε * (1 − λ)       , if r <= 2**(1/6) * σ
!            λ * u_LJ(r)                 , if r >  2**(1/6) * σ
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!
!   5. HPS_cation_pi (HPS-cation-pi-i and -ii): HPS_piecewise + LJ for specific π interactions
!      - U_total(r)        = U_HPS_piecewise(r) + U_extra_LJ(r) for cation-π pairs
!      - U_Electrostatic   = (q1 * q2 / r) * exp(−κ * r)
!
!   6. Pedone: Empirical potential for oxides/silicates
!      - U(r) = A * exp(−r / ρ) − C / r**6 + q1 * q2 / r
!
!   7. Tersoff: Bond-order potential for covalent systems
!      - U = 0.5 * ∑_i ∑_j≠i f_C(r_ij) * [ f_R(r_ij) + b_ij * f_A(r_ij) ]
!        where:
!          f_R(r) = A * exp(−λ1 * r)
!          f_A(r) = −B * exp(−λ2 * r)
!          b_ij   = bond order term (depends on angles and neighbors)
!
!   This module loads the appropriate force field modules that define:
!     - Intermolecular energy and pressure calculation routines
!     - Intramolecular CBMC/FECBMC Rosenbluth energy routines
!     - Trial move biasing (quick ΔE evaluators)
!
! Author:
!   Aliasghar Sepehri
!===================================================================
module EnergyPressurePointers

  !------------------------------------------------------
  ! Energy subroutine interfaces (used in MC trial evaluation)
  use E_Interface_LJ_Q
  use E_Interface_Mpipi
  use E_Interface_HPS_single
  use E_Interface_HPS_piecewise
  use E_Interface_HPS_cation_pi
  use E_Interface_Pedone
  use E_Interface_Tersoff

  !------------------------------------------------------
  ! Rosenbluth sampling routines (used in CBMC and AVBMC regrowth)
  use Rosenbluth_Functions_LJ_Q
  use Rosenbluth_Functions_Mpipi
  use Rosenbluth_Functions_HPS_single
  use Rosenbluth_Functions_HPS_piecewise
  use Rosenbluth_Functions_HPS_cation_pi
  use Rosenbluth_Functions_Pedone
  use Rosenbluth_Functions_Tersoff

  !------------------------------------------------------
  ! Pressure calculation subroutines (Virial-based)
  use Pressure_LJ_Electro
  use Pressure_Mpipi
  use Pressure_HPS_single
  use Pressure_HPS_piecewise
  use Pressure_HPS_cation_pi

  !------------------------------------------------------
  ! Quick neighbor energy evaluation (for AVBMC or pre-screening)
  use InterEnergy_LJ_Electro,       only: QuickNei_ECalc_Inter_LJ_Q
  use InterEnergy_Mpipi,            only: QuickNei_ECalc_Inter_Mpipi
  use InterEnergy_HPS_single,       only: QuickNei_ECalc_Inter_HPS_single
  use InterEnergy_HPS_piecewise,    only: QuickNei_ECalc_Inter_HPS_piecewise
  use InterEnergy_HPS_cation_pi,    only: QuickNei_ECalc_Inter_HPS_cation_pi
  use InterEnergy_Pedone,           only: QuickNei_ECalc_Inter_Pedone
  use InterEnergy_Tersoff,          only: QuickNei_ECalc_Inter_Tersoff
  !------------------------------------------------------
  ! Interface: DetailedEnergyInterface
  ! Purpose:
  !   Unified interface for computing the total energy of the system
  !   using the selected force field. Typically called at the beginning
  !   or during a Monte Carlo move.
  !
  !   Implementations are force-field specific and must be declared
  !   in modules like E_Interface_Mpipi, E_Interface_HPS_piecewise, etc.
  !
  ! Parameters:
  !   E_T      : [inout] Real(dp)
  !             - Output: total energy of the system or trial configuration
  !
  !   rejMove  : [inout] Logical
  !             - Initially .false., set to .true. if the configuration
  !               is physically invalid (e.g., atomic overlap detected)
  !
  !   This flag is especially important when reading input configurations,
  !   where overlapping atoms must be rejected before simulation begins.
  !
  interface
    subroutine DetailedEnergyInterface(E_T, rejMove)
      use VarPrecision
      implicit none
      real(dp), intent(inout) :: E_T         ! Total energy to be computed
      logical , intent(inout) :: rejMove     ! Set true if overlap or invalid geometry
    end subroutine DetailedEnergyInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: ShiftEnergyInterface
  ! Purpose:
  !   Computes the energy difference (∆E) resulting from a trial move
  !   (e.g., molecule translation, rotation, CBMC, FECBMC). It handles
  !   both intra- and inter-molecular interactions and is force-field specific.
  !
  !   This interface is used inside the core MC engine to evaluate 
  !   acceptance criteria via the Metropolis rule.
  !
  ! Parameters:
  !
  !   E_Inter     : [out]   Real(dp)
  !                 Intermolecular energy of the trial configuration.
  !
  !   E_Intra     : [out]   Real(dp)
  !                 Intra-molecular energy of the displaced molecule.
  !
  !   disp(:)     : [in]    Type(Displacement)
  !                 Array of atoms that have been displaced in the move.
  !
  !   PairList(:) : [inout] Real(dp)
  !                 Stores interaction values between the moved molecule
  !                 and all other molecules in the cluster:
  !                   - For energy-based cluster criteria: interaction energy.
  !                   - For distance-based: distance between first atoms.
  !                   - For minimum distance-based: minimum atomic distance.
  !
  !   dETable(:)  : [inout] Real(dp)
  !                 Stores the change in inter-molecular interaction energy
  !                 between the moved molecule and each cluster molecule.
  !
  !   useIntra(1:4) : [in]  Logical
  !                   Flags to indicate whether to calculate:
  !                   (1) bond, (2) angle, (3) torsion, (4) nonbonded intra terms.
  !
  !   useInter    : [in]    Logical, optional
  !                 If present and .false., skips intermolecular energy calculation.
  !
  !   NeighborDetailsNew(:) : [inout] Type(NeighborDetails), optional
  !                 Updated neighbor tracking information for clustering algorithms.
  !
  !   rejMove     : [inout] Logical
  !                 Set to .true. if the move results in overlap or another invalid configuration.
  interface
    subroutine ShiftEnergyInterface(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove, &
                                     useInter, NeighborDetailsNew)
      use CoordinateTypes
      use VarPrecision
      implicit none

      ! Optional flags
      logical, intent(in), optional :: useInter
      type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)

      ! Mandatory inputs
      logical, intent(in) :: useIntra(1:4)
      type(Displacement), intent(in) :: disp(:)

      ! Outputs
      real(dp), intent(out) :: E_Inter, E_Intra

      ! In/out accumulators
      real(dp), intent(inout) :: PairList(:)
      real(dp), intent(inout) :: dETable(:)
      logical, intent(inout) :: rejMove

    end subroutine ShiftEnergyInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: SwapInEnergyInterface
  ! Purpose:
  !   Computes the intra- and inter-molecular energy contributions
  !   when inserting a new molecule into the cluster as part of a 
  !   grand canonical Swap-In move.
  !
  !   Used to evaluate acceptance probability based on:
  !     ΔE = E_trial - E_current = E_Intra + E_Inter
  !
  ! Parameters:
  !
  !   E_Inter     : [out] Real(dp)
  !                 Intermolecular energy of the inserted molecule with the cluster.
  !
  !   E_Intra     : [out] Real(dp)
  !                 Internal energy of the newly inserted molecule (bond, angle, etc.).
  !
  !   PairList(:) : [inout] Real(dp)
  !                 Stores interaction values between the inserted molecule
  !                 and all other molecules in the cluster:
  !                   - For energy-based cluster criteria: interaction energy.
  !                   - For distance-based: distance between first atoms.
  !                   - For minimum distance-based: minimum atomic distance.
  !
  !   dETable(:)  : [inout] Real(dp)
  !                 Stores the change in inter-molecular interaction energy
  !                 between the inserted molecule and each cluster molecule.
  !
  !   useInter    : [in] Logical (optional)
  !                 If present and .false., skip intermolecular calculation.
  !
  !   NeighborDetailsNew(:) : [inout] Type(NeighborDetails) (optional)
  !                 Optional structure to record/update new neighbor relationships.
  !
  !   rejMove     : [out] Logical
  !                 Set to .true. if insertion results in invalid configuration (e.g., overlaps).
  !
  interface
    subroutine SwapInEnergyInterface(E_Inter, E_Intra, PairList, dETable, rejMove, &
                                     useInter, NeighborDetailsNew)
      use CoordinateTypes
      use VarPrecision
      implicit none

      ! Outputs
      real(dp), intent(out) :: E_Inter, E_Intra
      logical, intent(out) :: rejMove

      ! Optional inputs
      logical, intent(in), optional :: useInter
      type(NeighborDetails), intent(inout), optional :: NeighborDetailsNew(:)

      ! Energy accumulators
      real(dp), intent(inout) :: PairList(:), dETable(:)

    end subroutine SwapInEnergyInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: SwapOutEnergyInterface
  ! Purpose:
  !   Computes the intra- and inter-molecular energy of a molecule that
  !   is being removed (deleted) from the cluster in a grand canonical
  !   Swap-Out Monte Carlo move.
  !
  !   This energy is needed to compute:
  !     ΔE = E_removed = E_Intra + E_Inter
  !
  ! Parameters:
  !
  !   E_Inter     : [out] Real(dp)
  !                 Intermolecular energy of the molecule being deleted.
  !
  !   E_Intra     : [out] Real(dp)
  !                 Internal energy (bonded/nonbonded) of the molecule.
  !
  !   nType       : [in] Integer
  !                 Molecule type of the one being removed.
  !
  !   nMol        : [in] Integer
  !                 Molecule index (within type) of the one being removed.
  !
  !   dETable(:)  : [inout] Real(dp)
  !                 Stores the change in inter-molecular interaction energy
  !                 between the deleted molecule and each cluster molecule.
  !
  !   useInter    : [in] Logical (optional)
  !                 If present and .false., inter-molecular interactions are skipped.
  !
  interface
    subroutine SwapOutEnergyInterface(E_Inter, E_Intra, nType, nMol, dETable, useInter)
      use VarPrecision
      implicit none

      ! Output: energy of molecule to be removed
      real(dp), intent(out) :: E_Inter, E_Intra

      ! Molecule identity
      integer, intent(in) :: nType, nMol

      ! Energy difference array
      real(dp), intent(inout) :: dETable(:)

      ! Optional: skip inter-molecular interactions
      logical, intent(in), optional :: useInter
    end subroutine SwapOutEnergyInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: DetailedPressureInterface
  ! Purpose:
  !   Generic interface to compute the pressure of the current cluster
  !   configuration. The implementation is force-field specific and should
  !   apply the appropriate pressure formula (e.g., Virial equation).
  !
  !   This is used when `calcPressure = .true.` in SimParameters and
  !   pressure is reported during or after the simulation.
  !
  ! Parameters:
  !
  !   P_T : [out] Real(dp)
  !         Instantaneous pressure calculated from the current configuration.
  !
  interface
    subroutine DetailedPressureInterface(P_T)
      use VarPrecision
      implicit none
      real(dp), intent(out) :: P_T
    end subroutine DetailedPressureInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: ShiftPressureInterface
  ! Purpose:
  !   Computes the pressure of a trial configuration after a displacement
  !   move, without modifying the actual configuration. This is useful for
  !   pressure-based biasing or diagnostics during MC sampling.
  !
  ! Parameters:
  !
  !   disp(:)    : [in] Type(Displacement)
  !                List of displaced atoms or molecules for the trial move.
  !
  !   P_Trial    : [out] Real(dp)
  !                Pressure of the trial configuration after displacement.
  !
  interface
    subroutine ShiftPressureInterface(P_Trial, disp)
      use CoordinateTypes
      use VarPrecision
      implicit none

      ! Trial displacement information
      type(Displacement), intent(in) :: disp(:)

      ! Output: pressure of the trial configuration
      real(dp), intent(out) :: P_Trial
    end subroutine ShiftPressureInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: SwapInPressureInterface
  ! Purpose:
  !   Computes the pressure of a trial configuration after a molecule
  !   has been inserted (Swap-In move). Used in grand canonical MC
  !   simulations when pressure-sensitive biasing or logging is enabled.
  !
  ! Parameters:
  !
  !   P_Trial : [out] Real(dp)
  !             Pressure of the configuration with the inserted molecule.
  !
  interface
    subroutine SwapInPressureInterface(P_Trial)
      use VarPrecision
      implicit none

      ! Output: pressure after insertion
      real(dp), intent(out) :: P_Trial
    end subroutine SwapInPressureInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: SwapOutPressureInterface
  ! Purpose:
  !   Computes the pressure of the trial configuration after deleting
  !   a specific molecule during a Swap-Out (grand canonical deletion) move.
  !
  !   This is used in simulations where pressure feedback is monitored or 
  !   biasing depends on pressure after particle removal.
  !
  ! Parameters:
  !
  !   iType   : [in] Integer
  !             Molecule type of the one being removed.
  !
  !   iMol    : [in] Integer
  !             Index of the molecule within that type to be removed.
  !
  !   P_Trial : [out] Real(dp)
  !             Pressure after removing the specified molecule.
  !
  interface
    subroutine SwapOutPressureInterface(iType, iMol, P_Trial)
      use VarPrecision
      implicit none

      integer, intent(in) :: iType, iMol
      real(dp), intent(out) :: P_Trial
    end subroutine SwapOutPressureInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: RosenMolNewInterface
  ! Purpose:
  !   Performs Rosenbluth sampling by generating an entire molecule
  !   in one step (i.e., without sequential growth). This is typically
  !   used for small molecules or rigid bodies.
  !
  !   Called in CBMC or AVBMC insertion/regrowth steps.
  !
  ! Parameters:
  !
  !   nRosen   : [in] Integer
  !              Number of Rosenbluth trials (independent full configurations)
  !
  !   nType    : [in] Integer
  !              Molecule type to be generated
  !
  !   included(:) : [in] Logical
  !       Boolean mask over molecules in the cluster:
  !         - `.true.` means include this molecule when computing interactions
  !         - Only one trial (new or old) is evaluated at a time,
  !           so this array is reset per trial
  !
  !   E_Trial  : [out] Real(dp)
  !              Energy of the selected configuration (after weighted sampling)
  !
  !   overlap  : [out] Logical
  !              Set to .true. if overlap was detected in all trial configurations
  !
  interface
    subroutine RosenMolNewInterface(nRosen, nType, included, E_Trial, overlap)
      use VarPrecision
      implicit none

      ! Inputs
      integer, intent(in) :: nRosen, nType
      logical, intent(in) :: included(:)

      ! Outputs
      real(dp), intent(out) :: E_Trial
      logical, intent(out) :: overlap
    end subroutine RosenMolNewInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: RosenMolOldInterface
  ! Purpose:
  !   Computes the energy of the *existing (old)* molecular configuration
  !   during Rosenbluth sampling. This is used to calculate the acceptance
  !   probability in CBMC/AVBMC moves.
  !
  ! Parameters:
  !
  !   mol_x(:), mol_y(:), mol_z(:) : [in] Real(dp)
  !       Coordinates of atoms in the molecule to be regrown (current configuration)
  !
  !   nType     : [in] Integer
  !       Molecule type being regrown
  !
  !   included(:) : [in] Logical
  !       Boolean mask over molecules in the cluster:
  !         - `.true.` means include this molecule when computing interactions
  !         - Only one trial (new or old) is evaluated at a time,
  !           so this array is reset per trial
  !
  !   E_Trial   : [out] Real(dp)
  !       Energy of the current configuration computed using selected neighbors
  !
  interface
    subroutine RosenMolOldInterface(mol_x, mol_y, mol_z, nType, included, E_Trial)
      use VarPrecision
      implicit none

      integer, intent(in) :: nType
      logical, intent(in) :: included(:)            ! Neighbor molecule inclusion mask
      real(dp), intent(in) :: mol_x(:), mol_y(:), mol_z(:)
      real(dp), intent(out) :: E_Trial
    end subroutine RosenMolOldInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: RosenAtomNewInterface
  ! Purpose:
  !   Computes the interaction energy of a single atom trial position
  !   during Rosenbluth sampling in Configurational Bias Monte Carlo (CBMC)
  !   for **flexible small molecules**.
  !
  !   Atoms are grown sequentially (one at a time), 
  !
  ! Parameters:
  !
  !   nType     : [in] Integer
  !       Type of the molecule being regrown
  !
  !   nAtom     : [in] Integer
  !       Index of the atom currently being grown
  !
  !   trialPos  : [in] SimpleAtomCoords
  !       Coordinates of the trial atom position
  !
  !   included(:) : [in] Logical
  !       Mask indicating which other molecules in the system should be
  !       included in the energy evaluation for this trial.
  !       - Typically excludes the molecule currently being grown.
  !
  !   E_Trial   : [out] Real(dp)
  !       Computed interaction energy of this trial atom position
  !
  !   overlap   : [inout] Logical
  !       Set to `.true.` if this trial atom overlaps with a neighbor
  !
  interface
    subroutine RosenAtomNewInterface(nType, nAtom, trialPos, included, E_Trial, overlap)
      use CoordinateTypes
      use VarPrecision
      implicit none

      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      logical, intent(in) :: included(:)
      real(dp), intent(out) :: E_Trial
      logical, intent(inout) :: overlap
    end subroutine RosenAtomNewInterface
  end interface 
!======================================================
  !------------------------------------------------------
  ! Interface: RosenAtomOldInterface
  ! Purpose:
  !   Computes the interaction energy of the *existing atom configuration*
  !   for a specific atom within a flexible small molecule during CBMC.
  !
  !   This is used in computing the Rosenbluth weight for the "old" state
  !   during Metropolis acceptance of a regrowth move.
  !
  ! Parameters:
  !
  !   nType   : [in] Integer
  !       Type of the molecule being regrown
  !
  !   nMol    : [in] Integer
  !       Index of the molecule being regrown
  !
  !   nAtom   : [in] Integer
  !       Index of the atom currently being evaluated
  !
  !   included(:) : [in] Logical
  !       Mask specifying which other molecules in the system should be
  !       considered when computing the interaction energy.
  !       - Set `.false.` for distant or non-interacting molecules
  !
  !   E_Trial : [out] Real(dp)
  !       Computed energy of the current (old) atom position
  !
  interface
    subroutine RosenAtomOldInterface(nType, nMol, nAtom, included, E_Trial)
      use VarPrecision
      implicit none

      integer, intent(in) :: nType, nMol, nAtom
      logical, intent(in) :: included(:)
      real(dp), intent(out) :: E_Trial
    end subroutine RosenAtomOldInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: RosenAtomIntraNewInterface
  ! Purpose:
  !   Computes the nonbonded intramolecular interaction energy between
  !   a *trial atom position* and all **already regrown atoms** in the
  !   same molecule during CBMC regrowth.
  !
  !   This is necessary for flexible molecules where internal clashes
  !   or overlaps can arise between nonbonded atoms (e.g., via 1-4, 1-5
  !   interactions) that are not constrained by bonds/angles.
  !
  ! Parameters:
  !
  !   nType   : [in] Integer
  !       Molecule type being regrown
  !
  !   nAtom   : [in] Integer
  !       Index of the atom currently being placed
  !
  !   trialPos : [in] SimpleAtomCoords
  !       Coordinates of the trial position for this atom
  !
  !   regrown(:) : [in] Logical
  !       Boolean mask over atoms in the same molecule:
  !       - `.true.` for atoms already regrown and fixed in this CBMC step
  !       - Used to limit the internal interaction calculation
  !
  !   E_Trial : [out] Real(dp)
  !       Calculated nonbonded intramolecular energy contribution
  !       for this trial atom
  !
  !   overlap : [inout] Logical
  !       Initially set to `.false.`
  !       Set to `.true.` if the trial atom overlaps with a regrown atom
  !
  interface
    subroutine RosenAtomIntraNewInterface(nType, nAtom, trialPos, regrown, E_Trial, overlap)
      use CoordinateTypes
      use VarPrecision
      implicit none

      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      logical, intent(in) :: regrown(:)
      real(dp), intent(out) :: E_Trial
      logical, intent(inout) :: overlap
    end subroutine RosenAtomIntraNewInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: RosenAtomIntraOldInterface
  ! Purpose:
  !   Computes the nonbonded intramolecular interaction energy between the
  !   **existing position of an atom** and **already regrown atoms** within the
  !   same molecule, during Rosenbluth weight computation for CBMC.
  !
  !   Used to calculate the intramolecular contribution to the Rosenbluth
  !   factor of the "old" configuration.
  !
  ! Parameters:
  !
  !   nType   : [in] Integer
  !       Molecule type of the molecule being regrown
  !
  !   nMol    : [in] Integer
  !       Molecule index within the simulation
  !
  !   nAtom   : [in] Integer
  !       Atom index being evaluated
  !
  !   trialPos : [in] SimpleAtomCoords
  !       Coordinates of the atom in the *old* (existing) configuration
  !
  !   regrown(:) : [in] Logical
  !       Mask of atoms already regrown in the molecule
  !
  !   E_Trial : [out] Real(dp)
  !       Computed nonbonded intramolecular energy
  !       between the selected atom and previously regrown atoms
  !
  interface
    subroutine RosenAtomIntraOldInterface(nType, nMol, nAtom, trialPos, regrown, E_Trial)
      use CoordinateTypes
      use VarPrecision
      implicit none

      integer, intent(in) :: nType, nAtom, nMol
      type(SimpleAtomCoords), intent(in) :: trialPos
      logical, intent(in) :: regrown(:)
      real(dp), intent(out) :: E_Trial
    end subroutine RosenAtomIntraOldInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: RosenAtomIntraGasOldInterface
  ! Purpose:
  !   Computes the intramolecular (nonbonded) energy of an atom in a molecule
  !   located in the dilute or gas phase, during 
  !   insertion proposals in CBMC/AVBMC.
  !
  !   Used to ensure that intramolecular interactions are properly accounted for
  !   when a molecule is inserted from or deleted into the dilute reservoir.
  !
  ! Parameters:
  !
  !   nType   : [in] Integer
  !       Molecule type being inserted/deleted
  !
  !   nAtom   : [in] Integer
  !       Index of the atom being evaluated
  !
  !   trialPos : [in] SimpleAtomCoords
  !       Position of the atom in the gas-phase configuration
  !
  !   regrown(:) : [in] Logical
  !       Mask of atoms already placed in the molecule (used for computing
  !       internal energy contributions with the growing atom)
  !
  !   E_Trial : [out] Real(dp)
  !       Calculated nonbonded intramolecular energy for this atom
  !
  interface
    subroutine RosenAtomIntraGasOldInterface(nType, nAtom, trialPos, regrown, E_Trial)
      use CoordinateTypes
      use VarPrecision
      implicit none

      integer, intent(in) :: nType, nAtom
      type(SimpleAtomCoords), intent(in) :: trialPos
      logical, intent(in) :: regrown(:)
      real(dp), intent(out) :: E_Trial
    end subroutine RosenAtomIntraGasOldInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: RosenLLAtomNewInterface
  ! Purpose:
  !   Computes the interaction energy of a *trial atom position* with a selected
  !   set of atoms (typically neighbors), during regrowth of Long Linear (LL)
  !   biomolecules using CBMC or AVBMC.
  !
  !   This method is optimized for performance by only evaluating interactions
  !   against atoms listed in `InterAtms`, rather than checking all atoms in
  !   the cluster or full molecule.
  !
  ! Parameters:
  !
  !   nType   : [in] Integer
  !       Type of the molecule being regrown
  !
  !   nAtom   : [in] Integer
  !       Index of the atom being added/regrown
  !
  !   trialPos : [in] SimpleAtomCoords
  !       Coordinates of the trial position of the atom
  !
  !   nInterAtms : [in] Integer
  !       Number of atoms to evaluate interactions with
  !
  !   InterAtms(3, nInterAtms) : [in] Integer
  !       A 2D array containing identifiers of atoms to interact with.
  !       Structured as:
  !         - InterAtms(1,i): Molecule type
  !         - InterAtms(2,i): Molecule index
  !         - InterAtms(3,i): Atom index in that molecule
  !
  !   E_Trial : [out] Real(dp)
  !       Total interaction energy between the trial atom and selected neighbors
  !
  !   overlap : [inout] Logical
  !       Set to `.true.` if any overlap is detected with the neighbors
  !
  interface
    subroutine RosenLLAtomNewInterface(nType, nAtom, trialPos, nInterAtms, InterAtms, E_Trial, overlap)
      use CoordinateTypes
      use VarPrecision
      implicit none

      integer, intent(in) :: nType, nAtom, nInterAtms
      type(SimpleAtomCoords), intent(in) :: trialPos
      integer, intent(in) :: InterAtms(:,:)
      real(dp), intent(out) :: E_Trial
      logical, intent(inout) :: overlap
    end subroutine RosenLLAtomNewInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: RosenLLAtomOldInterface
  ! Purpose:
  !   Computes the interaction energy of an existing atom in a Long Linear (LL)
  !   biomolecule with a selected set of atoms. Used in:
  !     - Rosenbluth weight computation during deletion or regrowth
  !     - Evaluating contribution of the old atom before replacement
  !
  ! Parameters:
  !
  !   nType   : [in] Integer
  !       Molecule type of the atom
  !
  !   nMol    : [in] Integer
  !       Molecule index (i.e., which instance of the molecule type)
  !
  !   nAtom   : [in] Integer
  !       Index of the atom in the molecule
  !
  !   nInterAtms : [in] Integer
  !       Number of atoms to evaluate interactions with
  !
  !   InterAtms(3, nInterAtms) : [in] Integer
  !       2D array listing atoms to consider for interaction.
  !       Structured as:
  !         - InterAtms(1,i): Molecule type
  !         - InterAtms(2,i): Molecule index
  !         - InterAtms(3,i): Atom index in that molecule
  !
  !   E_Trial : [out] Real(dp)
  !       Total interaction energy of the specified atom with the listed neighbors
  !
  interface
    subroutine RosenLLAtomOldInterface(nType, nMol, nAtom, nInterAtms, InterAtms, E_Trial)
      use VarPrecision
      implicit none

      integer, intent(in) :: nType, nMol, nAtom, nInterAtms
      integer, intent(in) :: InterAtms(:,:)
      real(dp), intent(out) :: E_Trial
    end subroutine RosenLLAtomOldInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: RosenLLAtomIntraNewInterface
  !
  ! Purpose:
  !   Computes the **intra-molecular** (nonbonded) interaction energy of a 
  !   trial atom position with a specified list of atoms from the **same molecule**
  !   during CBMC/AVBMC regrowth of long linear (LL) biomolecules.
  !
  !   Designed for flexible chains, where bonded atoms are excluded and 
  !   only relevant intra-molecular pairs are evaluated.
  !
  ! Parameters:
  !
  !   nType : [in] Integer
  !       Molecule type
  !
  !   nAtom : [in] Integer
  !       Index of the atom being regrown within the molecule
  !
  !   trialPos : [in] SimpleAtomCoords
  !       Coordinates of the trial atom
  !
  !   nIntraAtms : [in] Integer
  !       Number of intra-molecular atoms to check interactions with
  !
  !   IntraAtms(nIntraAtms) : [in] Integer
  !       List of atom indices within the same molecule to check for interactions
  !
  !   E_Trial : [out] Real(dp)
  !       Calculated intra-molecular energy of the trial atom
  !
  !   overlap : [inout] Logical
  !       Set to `.true.` if the trial atom overlaps with any intra-molecular atom
  !
  interface
    subroutine RosenLLAtomIntraNewInterface(nType, nAtom, trialPos, nIntraAtms, IntraAtms, E_Trial, overlap)
      use CoordinateTypes
      use VarPrecision
      implicit none

      integer, intent(in) :: nType, nAtom, nIntraAtms
      integer, intent(in) :: IntraAtms(:)
      type(SimpleAtomCoords), intent(in) :: trialPos
      real(dp), intent(out) :: E_Trial
      logical, intent(inout) :: overlap
    end subroutine RosenLLAtomIntraNewInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: RosenLLAtomIntraOldInterface
  !
  ! Purpose:
  !   Computes the **intra-molecular** (nonbonded) interaction energy of an
  !   existing atom (already in the molecule) with a list of other atoms 
  !   within the **same** LL biomolecule.
  !
  !   Used during CBMC/AVBMC regrowth to evaluate the contribution of the
  !   atom to the original Rosenbluth weight.
  !
  ! Parameters:
  !
  !   nType : [in] Integer
  !       Molecule type of the atom
  !
  !   nMol : [in] Integer
  !       Index of the molecule
  !
  !   nAtom : [in] Integer
  !       Index of the atom within the molecule
  !
  !   trialPos : [in] SimpleAtomCoords
  !       Position to use for the atom in energy calculations 
  !       (can be current or trial during hybrid reuse)
  !
  !   nIntraAtms : [in] Integer
  !       Number of atoms to evaluate interactions with
  !
  !   IntraAtms(nIntraAtms) : [in] Integer
  !       List of intra-molecular atom indices within the same molecule
  !
  !   E_Trial : [out] Real(dp)
  !       Computed intra-molecular energy contribution from this atom
  !
  interface
    subroutine RosenLLAtomIntraOldInterface(nType, nMol, nAtom, trialPos, nIntraAtms, IntraAtms, E_Trial)
      use CoordinateTypes
      use VarPrecision
      implicit none

      integer, intent(in) :: nType, nMol, nAtom, nIntraAtms
      integer, intent(in) :: IntraAtms(:)
      type(SimpleAtomCoords), intent(in) :: trialPos
      real(dp), intent(out) :: E_Trial
    end subroutine RosenLLAtomIntraOldInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: RosenLLAtomIntraGasOldInterface
  !
  ! Purpose:
  !   Calculates the **intra-molecular** (nonbonded) energy of an atom in 
  !   a **long linear (LL)** biomolecule in its **original gas-phase state** 
  !   *before* a SwapIn (insertion) Monte Carlo move.
  !
  !   This subroutine is used to compute the old-state energy contribution 
  !   from the atom’s intra-molecular (same molecule) interactions while it 
  !   is still in the dilute (gas) phase.
  !
  !   It is used when computing the Rosenbluth factor for CBMC/AVBMC insertion.
  !
  ! Parameters:
  !
  !   nType : [in] Integer  
  !       Type of the molecule that the atom belongs to
  !
  !   nAtom : [in] Integer  
  !       Atom index within the molecule
  !
  !   trialPos : [in] SimpleAtomCoords  
  !       Current position of the atom in the gas phase
  !
  !   nIntraAtms : [in] Integer  
  !       Number of atoms in the same molecule used for intra-energy evaluation
  !
  !   IntraAtms(nIntraAtms) : [in] Integer  
  !       List of intra-molecular atom indices for energy calculation
  !
  !   E_Trial : [out] Real(dp)  
  !       Computed intra-molecular energy of the atom in its old gas-phase state
  !
  interface
    subroutine RosenLLAtomIntraGasOldInterface(nType, nAtom, trialPos, nIntraAtms, IntraAtms, E_Trial)
      use CoordinateTypes
      use VarPrecision
      implicit none

      integer, intent(in) :: nType, nAtom, nIntraAtms
      integer, intent(in) :: IntraAtms(:)
      type(SimpleAtomCoords), intent(in) :: trialPos
      real(dp), intent(out) :: E_Trial
    end subroutine RosenLLAtomIntraGasOldInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: FindInteractingAtomsInterface
  !
  ! Purpose:
  !   Identifies atoms that are within a specified interaction region defined
  !   as the **cutoff distance plus a buffer radius `r0`** from a given 
  !   reference atom. This is typically used in CBMC/AVBMC regrowth steps 
  !   of long linear (LL) biomolecules to locate potential interaction partners.
  !
  !   Returns a list of interacting atoms as a 2D array with each entry 
  !   identifying molecule type, molecule index, and atom index.
  !
  ! Parameters:
  !
  !   RefAtm : [in] SimpleAtomCoords  
  !       Cartesian coordinates (x, y, z) of the reference atom
  !
  !   r0 : [in] Real(dp)  
  !       Buffer distance to be added to the force field cutoff when 
  !       determining which atoms to include as potential interactors
  !
  !   iType : [in] Integer(kind=atomIntType)  
  !       Type of the reference molecule (used to exclude self-interactions)
  !
  !   nIndx : [in] Integer  
  !       Index of the reference molecule (used to exclude self-interactions)
  !
  !   nInterAtoms : [inout] Integer  
  !       On input: ignored or set to 0  
  !       On output: number of interacting atoms found
  !
  !   InterAtoms(:,:) : [out] Integer  
  !       List of interacting atoms in the format:
  !         - InterAtoms(1,i): Molecule type
  !         - InterAtoms(2,i): Molecule index
  !         - InterAtoms(3,i): Atom index in that molecule
  !
  interface
    subroutine FindInteractingAtomsInterface(RefAtm, r0, iType, nIndx, nInterAtoms, InterAtoms)
      use CoordinateTypes
      use VarPrecision
      implicit none

      type(SimpleAtomCoords), intent(in) :: RefAtm
      real(dp), intent(in) :: r0
      integer(kind=atomIntType), intent(in) :: iType
      integer, intent(in) :: nIndx
      integer, intent(inout) :: nInterAtoms
      integer, intent(out) :: InterAtoms(:,:)
    end subroutine FindInteractingAtomsInterface
  end interface
!======================================================
  !------------------------------------------------------
  ! Interface: QuickInterface
  !
  ! Purpose:
  !   This subroutine is used during particle **insertion moves** 
  !   (SwapIn) to quickly determine whether the inserted molecule
  !   is connected to the existing cluster based on an 
  !   **energy-based cluster criterion**.
  !
  !   Specifically, it checks whether the newly inserted molecule
  !   has a **nonzero interaction energy** with **at least one**
  !   molecule already in the cluster. If not, the move is rejected.
  !
  ! Parameters:
  !
  !   jType : [in] Integer  
  !       Type of the molecule being inserted
  !
  !   jMol : [in] Integer  
  !       Index of the molecule being inserted
  !
  !   rejMove : [out] Logical  
  !       Set to `.true.` if the molecule is **not connected** to the cluster
  !       (i.e., no interaction energy with any existing molecule)
  !
  interface
    subroutine QuickInterface(jType, jMol, rejMove)
      implicit none
      integer, intent(in) :: jType, jMol
      logical, intent(out) :: rejMove
    end subroutine QuickInterface
  end interface
!======================================================

!============================== Force Field Interface Procedure Pointers ==============================
! These procedure pointers allow the Monte Carlo simulation to dynamically link to the appropriate
! force field-dependent subroutines for energy, pressure, and Rosenbluth factor calculations.
! Each pointer is assigned at runtime based on user-defined force field in the simulation input.
!===============================================================================================

!========================= Energy Calculation Interfaces =========================
! Total energy evaluation for the current configuration (used for initial setup or overlap check)
procedure(DetailedEnergyInterface), pointer :: Detailed_ECalc => NULL()

! Energy change due to a molecular displacement (e.g., translation, rotation)
procedure(ShiftEnergyInterface), pointer  :: Shift_ECalc => NULL()

! Energy evaluation during molecule insertion (SwapIn move)
procedure(SwapInEnergyInterface), pointer :: SwapIn_ECalc => NULL()

! Energy evaluation during molecule deletion (SwapOut move)
procedure(SwapOutEnergyInterface), pointer :: SwapOut_ECalc => NULL()

!========================= Pressure Calculation Interfaces ========================
! Total pressure evaluation from all molecules (Virial contribution)
procedure(DetailedPressureInterface), pointer :: Detailed_PCalc => NULL()

! Pressure calculation during displacement move (trial configuration)
procedure(ShiftPressureInterface), pointer  :: Shift_PCalc => NULL()

! Pressure contribution from insertion (SwapIn)
procedure(SwapInPressureInterface), pointer :: SwapIn_PCalc => NULL()

! Pressure contribution removed during deletion (SwapOut)
procedure(SwapOutPressureInterface), pointer :: SwapOut_PCalc => NULL()

!========================= Rosenbluth Interfaces – Small Molecules =================
! Whole-molecule regrowth (new configuration): energy of inserted state
procedure(RosenMolNewInterface), pointer :: Rosen_Mol_New => NULL()

! Whole-molecule energy for old configuration (used for Rosenbluth factor calculation)
procedure(RosenMolOldInterface), pointer :: Rosen_Mol_Old => NULL()

!========================= Rosenbluth Interfaces – Atom-Based (Flexible) ===========
! Intermolecular energy of trial atom (new position)
procedure(RosenAtomNewInterface), pointer :: Rosen_Atom_New => NULL()

! Intermolecular energy of current atom position
procedure(RosenAtomOldInterface), pointer :: Rosen_Atom_Old => NULL()

! Intramolecular energy for a newly placed trial atom
procedure(RosenAtomIntraNewInterface), pointer :: Rosen_Atom_Intra_New => NULL()

! Intramolecular energy of atom in current configuration
procedure(RosenAtomIntraOldInterface), pointer :: Rosen_Atom_Intra_Old => NULL()

! Intramolecular energy in gas phase before insertion (old state)
procedure(RosenAtomIntraGasOldInterface), pointer :: Rosen_Atom_Intra_Gas_Old => NULL()

!========================= Rosenbluth Interfaces – Long Linear (LL) Molecules =======
! Intermolecular energy of trial atom for LL molecule (new position)
procedure(RosenLLAtomNewInterface), pointer :: RosenLL_Atom_New => NULL()

! Intermolecular energy of existing LL atom
procedure(RosenLLAtomOldInterface), pointer :: RosenLL_Atom_Old => NULL()

! Intramolecular energy for LL trial atom
procedure(RosenLLAtomIntraNewInterface), pointer :: RosenLL_Atom_Intra_New => NULL()

! Intramolecular energy for LL atom in original configuration
procedure(RosenLLAtomIntraOldInterface), pointer :: RosenLL_Atom_Intra_Old => NULL()

! Intramolecular energy in gas phase for LL molecule before insertion
procedure(RosenLLAtomIntraGasOldInterface), pointer :: RosenLL_Atom_Intra_Gas_Old => NULL()

!========================= Auxiliary and Cluster Criteria ==========================
! Find atoms within a certain distance of a reference atom (used in LL regrowth)
procedure(FindInteractingAtomsInterface), pointer :: Find_InterAtms => NULL()

! Quick rejection check based on energy cluster criterion (ensures connectivity)
procedure(QuickInterface), pointer :: Quick_Nei_ECalc => NULL()



end module EnergyPressurePointers
