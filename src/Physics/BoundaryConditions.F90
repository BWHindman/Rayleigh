!
!  Copyright (C) 2018 by the authors of the RAYLEIGH code.
!
!  This file is part of RAYLEIGH.
!
!  RAYLEIGH is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3, or (at your option)
!  any later version.
!
!  RAYLEIGH is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with RAYLEIGH; see the file LICENSE.  If not see
!  <http://www.gnu.org/licenses/>.
!

Module BoundaryConditions
    Use Math_Constants
    Use ProblemSize
    Use Fields
    Use Load_Balance, Only : mp_lm_values, l_lm_values, my_num_lm, &
                             m_lm_values, my_lm_min, my_nl_lm, my_nm_lm, &
                             my_lm_lval, my_lm_max, lm_count, lm_owner

    Use Generic_Input
    Use PDE_Coefficients
    Use Checkpointing, Only : Store_BC_Mask, Load_BC_Mask
    Use Controls, Only : full_restart

    Implicit None

    Logical :: Fix_Tvar_Top    = .False.
    Logical :: Fix_Tvar_Bottom = .False.
    Logical :: Fix_dTdr_Top    = .False.
    Logical :: Fix_dTdr_Bottom = .False.

    Logical :: Couple_Tvar_Top    = .False.
    Logical :: Couple_Tvar_Bottom = .False.
    Logical :: Couple_dTdr_Top    = .False.
    Logical :: Couple_dTdr_Bottom = .False.

    Logical :: Fix_chivar_a_Top(1:n_scalar_max)    = .False.
    Logical :: Fix_chivar_a_Bottom(1:n_scalar_max) = .False.
    Logical :: Fix_dchidr_a_Top(1:n_scalar_max)    = .False.
    Logical :: Fix_dchidr_a_Bottom(1:n_scalar_max) = .False.

    Logical :: Couple_chivar_a_Top(1:n_scalar_max)    = .False.
    Logical :: Couple_chivar_a_Bottom(1:n_scalar_max) = .False.
    Logical :: Couple_dchidr_a_Top(1:n_scalar_max)    = .False.
    Logical :: Couple_dchidr_a_Bottom(1:n_scalar_max) = .False.

    Logical :: Fix_chivar_p_Top(1:n_scalar_max)    = .False.
    Logical :: Fix_chivar_p_Bottom(1:n_scalar_max) = .False.
    Logical :: Fix_dchidr_p_Top(1:n_scalar_max)    = .False.
    Logical :: Fix_dchidr_p_Bottom(1:n_scalar_max) = .False.

    Logical :: Fix_divrt_top = .False.
    Logical :: Fix_divt_top = .False.
    Logical :: Fix_divrfc_top = .False.
    Logical :: Fix_divfc_top = .False.
    Logical :: Fix_poloidalfield_top = .False.
    Logical :: Fix_poloidalfield_bottom = .False.
    Logical :: Impose_Dipole_Field = .False.
    Logical :: Dipole_Field_Bottom = .False.
    Logical :: adjust_dTdr_Top = .True.

    Real*8  :: T_Bottom     = 0.0d0
    Real*8  :: T_Top        = 0.0d0
    Real*8  :: dTdr_Top     = 0.0d0
    Real*8  :: dTdr_Bottom  = 0.0d0

    Real*8  :: T_chi_a_coeff_top(1:n_scalar_max)     = 0.0d0
    Real*8  :: T_dchidr_a_coeff_top(1:n_scalar_max)  = 0.0d0
    Real*8  :: T_dTdr_coeff_top                    = 0.0d0
    Real*8  :: T_chi_a_coeff_bottom(1:n_scalar_max)     = 0.0d0
    Real*8  :: T_dchidr_a_coeff_bottom(1:n_scalar_max)  = 0.0d0
    Real*8  :: T_dTdr_coeff_bottom                    = 0.0d0
    Real*8  :: dTdr_chi_a_coeff_top(1:n_scalar_max)     = 0.0d0
    Real*8  :: dTdr_dchidr_a_coeff_top(1:n_scalar_max)  = 0.0d0
    Real*8  :: dTdr_T_coeff_top                       = 0.0d0
    Real*8  :: dTdr_chi_a_coeff_bottom(1:n_scalar_max)     = 0.0d0
    Real*8  :: dTdr_dchidr_a_coeff_bottom(1:n_scalar_max)  = 0.0d0
    Real*8  :: dTdr_T_coeff_bottom                       = 0.0d0

    Real*8  :: chi_a_Bottom(1:n_scalar_max)     = 0.0d0
    Real*8  :: chi_a_Top(1:n_scalar_max)        = 0.0d0
    Real*8  :: dchidr_a_Top(1:n_scalar_max)     = 0.0d0
    Real*8  :: dchidr_a_Bottom(1:n_scalar_max)  = 0.0d0

    Real*8  :: chi_a_chi_a_coeff_top(1:n_scalar_max,1:n_scalar_max)     = 0.0d0
    Real*8  :: chi_a_dchidr_a_coeff_top(1:n_scalar_max,1:n_scalar_max)  = 0.0d0
    Real*8  :: chi_a_T_coeff_top(1:n_scalar_max)                      = 0.0d0
    Real*8  :: chi_a_dTdr_coeff_top(1:n_scalar_max)                   = 0.0d0
    Real*8  :: chi_a_chi_a_coeff_bottom(1:n_scalar_max,1:n_scalar_max)     = 0.0d0
    Real*8  :: chi_a_dchidr_a_coeff_bottom(1:n_scalar_max,1:n_scalar_max)  = 0.0d0
    Real*8  :: chi_a_T_coeff_bottom(1:n_scalar_max)                      = 0.0d0
    Real*8  :: chi_a_dTdr_coeff_bottom(1:n_scalar_max)                   = 0.0d0
    Real*8  :: dchidr_a_chi_a_coeff_top(1:n_scalar_max,1:n_scalar_max)     = 0.0d0
    Real*8  :: dchidr_a_dchidr_a_coeff_top(1:n_scalar_max,1:n_scalar_max)  = 0.0d0
    Real*8  :: dchidr_a_T_coeff_top(1:n_scalar_max)                      = 0.0d0
    Real*8  :: dchidr_a_dTdr_coeff_top(1:n_scalar_max)                   = 0.0d0
    Real*8  :: dchidr_a_chi_a_coeff_bottom(1:n_scalar_max,1:n_scalar_max)     = 0.0d0
    Real*8  :: dchidr_a_dchidr_a_coeff_bottom(1:n_scalar_max,1:n_scalar_max)  = 0.0d0
    Real*8  :: dchidr_a_T_coeff_bottom(1:n_scalar_max)                      = 0.0d0
    Real*8  :: dchidr_a_dTdr_coeff_bottom(1:n_scalar_max)                   = 0.0d0

    Real*8  :: chi_p_Bottom(1:n_scalar_max)     = 1.0d0
    Real*8  :: chi_p_Top(1:n_scalar_max)        = 0.0d0
    Real*8  :: dchidr_p_Top(1:n_scalar_max)     = 0.0d0
    Real*8  :: dchidr_p_Bottom(1:n_scalar_max)  = 0.0d0

    Real*8  :: C10_bottom = 0.0d0
    Real*8  :: C10_top = 0.0d0
    Real*8  :: C11_bottom = 0.0d0
    Real*8  :: C11_top = 0.0d0
    Real*8  :: C1m1_bottom = 0.0d0
    Real*8  :: C1m1_top = 0.0d0
    Real*8  :: Br_bottom = 0.0d0
    Real*8  :: Dipole_Tilt_Degrees = 0.0d0

    Character*120 :: T_top_file         = '__nothing__'
    Character*120 :: T_bottom_file      = '__nothing__'
    Character*120 :: dTdr_top_file      = '__nothing__'
    Character*120 :: dTdr_bottom_file   = '__nothing__'
    Character*120 :: C_top_file         = '__nothing__'
    Character*120 :: C_bottom_file      = '__nothing__'

    Character*120 :: chi_a_top_file(1:n_scalar_max)       = '__nothing__'
    Character*120 :: chi_a_bottom_file(1:n_scalar_max)    = '__nothing__'
    Character*120 :: dchidr_a_top_file(1:n_scalar_max)    = '__nothing__'
    Character*120 :: dchidr_a_bottom_file(1:n_scalar_max) = '__nothing__'

    Character*120 :: chi_p_top_file(1:n_scalar_max)       = '__nothing__'
    Character*120 :: chi_p_bottom_file(1:n_scalar_max)    = '__nothing__'
    Character*120 :: dchidr_p_top_file(1:n_scalar_max)    = '__nothing__'
    Character*120 :: dchidr_p_bottom_file(1:n_scalar_max) = '__nothing__'

    Logical :: Strict_L_Conservation = .false.  ! Turn on to enforce angular momentum conservation 
                                                ! about x,y, and z-axes via integral constraint on Z(ell=1)
    Logical :: no_slip_boundaries = .false. ! Set to true to use no-slip boundaries.  Stree-free boundaries are the default.
    Logical :: stress_free_top = .true., stress_free_bottom = .true.
    Logical :: no_slip_top = .false., no_slip_bottom = .false.
    Logical :: use_checkpoint_bc_file = .false.

    Real*8, allocatable, dimension(:,:,:,:) :: bc_values  ! a 4-D array: (top/bottom, real/imag, my_num_lm, n_equations)

    Namelist /Boundary_Conditions_Namelist/ Fix_Tvar_Top, Fix_Tvar_Bottom, T_Bottom, T_Top, dTdr_top, dTdr_bottom, &
        fix_dtdr_bottom, fix_dtdr_top, fix_divrt_top, fix_divt_top, fix_divrfc_top, fix_divfc_top, &
        no_slip_boundaries, strict_L_Conservation, fix_poloidalfield_top, fix_poloidalfield_bottom, &
        C10_bottom, C10_top, C11_bottom, C11_top, C1m1_bottom, C1m1_top, Br_bottom, &
        dipole_tilt_degrees, impose_dipole_field, no_slip_top, no_slip_bottom, &
        stress_free_top, stress_free_bottom, T_top_file, T_bottom_file, dTdr_top_file, dTdr_bottom_file, &
        C_top_file, C_bottom_file, dipole_field_bottom, adjust_dTdr_top, &
        chi_a_top_file, chi_a_bottom_file, dchidr_a_top_file, dchidr_a_bottom_file, &
        fix_chivar_a_top, fix_chivar_a_bottom, chi_a_bottom, chi_a_top, &
        fix_dchidr_a_bottom, fix_dchidr_a_top, dchidr_a_top, dchidr_a_bottom, &
        chi_p_top_file, chi_p_bottom_file, dchidr_p_top_file, dchidr_p_bottom_file, &
        fix_chivar_p_top, fix_chivar_p_bottom, chi_p_bottom, chi_p_top, &
        fix_dchidr_p_bottom, fix_dchidr_p_top, dchidr_p_top, dchidr_p_bottom, &
        couple_Tvar_top, couple_Tvar_bottom, couple_dtdr_top, couple_dtdr_bottom, &
        T_chi_a_coeff_top, T_dchidr_a_coeff_top, T_dTdr_coeff_top, &
        T_chi_a_coeff_bottom, T_dchidr_a_coeff_bottom, T_dTdr_coeff_bottom, &
        dTdr_chi_a_coeff_top, dTdr_dchidr_a_coeff_top, dTdr_T_coeff_top, &
        dTdr_chi_a_coeff_bottom, dTdr_dchidr_a_coeff_bottom, dTdr_T_coeff_bottom, &
        couple_chivar_a_top, couple_chivar_a_bottom, couple_dchidr_a_top, couple_dchidr_a_bottom, &
        chi_a_chi_a_coeff_top, chi_a_dchidr_a_coeff_top, chi_a_T_coeff_top, chi_a_dTdr_coeff_top, &
        chi_a_chi_a_coeff_bottom, chi_a_dchidr_a_coeff_bottom, chi_a_T_coeff_bottom, chi_a_dTdr_coeff_bottom, &
        dchidr_a_chi_a_coeff_top, dchidr_a_dchidr_a_coeff_top, dchidr_a_T_coeff_top, dchidr_a_dTdr_coeff_top, &
        dchidr_a_chi_a_coeff_bottom, dchidr_a_dchidr_a_coeff_bottom, dchidr_a_T_coeff_bottom, dchidr_a_dTdr_coeff_bottom

Contains

    Subroutine Initialize_Boundary_Conditions()
        Implicit None
        Real*8 :: tilt_angle_radians,a,b
        Real*8 :: Ftop, rhotk_top, Fbottom
        Integer :: i
        Integer :: n_active_bcs

        n_active_bcs = count( (/ fix_tvar_top, fix_dtdr_top, couple_tvar_top, couple_dtdr_top /) )
        if (n_active_bcs > 1) then
           call stdout%print(" -- Error: Incompatible boundary conditions for tvar on top boundary.")
           call stdout%print("        ")      
           call pfi%exit(1)
        else if (n_active_bcs == 0) then
           call stdout%print(" -- Warning: no boundary conditions set for tvar on top boundary.")
           call stdout%print("        Defaulting to fix_tvar_top with T_top = 0.0.")
           call stdout%print("        ")      
           fix_tvar_top = .true.
           T_top = 0.0
        endif

        n_active_bcs = count( (/ fix_tvar_bottom, fix_dtdr_bottom, couple_tvar_bottom, couple_dtdr_bottom /) )
        if (n_active_bcs > 1) then
           call stdout%print(" -- Error: Incompatible boundary conditions for tvar on bottom boundary.")
           call stdout%print("        ")      
           call pfi%exit(1)
        else if (n_active_bcs == 0) then
           call stdout%print(" -- Warning: no boundary conditions set for tvar on bottom boundary.")
           call stdout%print("        Defaulting to fix_tvar_bottom with T_bottom = 1.0.")
           call stdout%print("        ")      
           fix_tvar_bottom = .true.
           T_bottom = 1.0
        endif

        do i = 1, n_active_scalars
          n_active_bcs = count( (/ fix_chivar_a_top(i), fix_dchidr_a_top(i), couple_chivar_a_top(i), couple_dchidr_a_top(i) /) )
          if (n_active_bcs > 1) then
             call stdout%print(" -- Error: Incompatible boundary conditions for chivar_a on top boundary.")
             call stdout%print("        ")      
             call pfi%exit(1)
          else if (n_active_bcs == 0) then
             call stdout%print(" -- Warning: no boundary conditions set for chivar_a on top boundary.")
             call stdout%print("        Defaulting to fix_chivar_a_top.")
             call stdout%print("        ")      
             fix_chivar_a_top(i) = .true.
          endif

          n_active_bcs = count( (/ fix_chivar_a_bottom(i), fix_dchidr_a_bottom(i), &
                                   couple_chivar_a_bottom(i), couple_dchidr_a_bottom(i) /) )
          if (n_active_bcs > 1) then
             call stdout%print(" -- Error: Incompatible boundary conditions for chivar_a on bottom boundary.")
             call stdout%print("        ")      
             call pfi%exit(1)
          else if (n_active_bcs == 0) then
             call stdout%print(" -- Warning: no boundary conditions set for chivar_a on bottom boundary.")
             call stdout%print("        Defaulting to fix_chivar_a_bottom.")
             call stdout%print("        ")      
             fix_chivar_a_bottom(i) = .true.
          endif
        end do

        do i = 1, n_passive_scalars
          n_active_bcs = count( (/ fix_chivar_p_top(i), fix_dchidr_p_top(i) /) )
          if (n_active_bcs > 1) then
             call stdout%print(" -- Error: Incompatible boundary conditions for chivar_p on top boundary.")
             call stdout%print("        ")      
             call pfi%exit(1)
          else if (n_active_bcs == 0) then
             call stdout%print(" -- Warning: no boundary conditions set for chivar_a on top boundary.")
             call stdout%print("        Defaulting to fix_chivar_p_top.")
             call stdout%print("        ")      
             fix_chivar_p_top(i) = .true.
          endif

          n_active_bcs = count( (/ fix_chivar_p_bottom(i), fix_dchidr_p_bottom(i) /) )
          if (n_active_bcs > 1) then
             call stdout%print(" -- Error: Incompatible boundary conditions for chivar_p on bottom boundary.")
             call stdout%print("        ")      
             call pfi%exit(1)
          else if (n_active_bcs == 0) then
             call stdout%print(" -- Warning: no boundary conditions set for chivar_a on bottom boundary.")
             call stdout%print("        Defaulting to fix_chivar_p_bottom.")
             call stdout%print("        ")      
             fix_chivar_p_bottom(i) = .true.
          endif
        end do

        If (no_slip_top .and. strict_L_Conservation) Then
            If (my_rank .eq. 0) Then
                Call stdout%print(" -- Warning:  Incompatible boundary conditions.")
                Call stdout%print("        Both strict_L_Conservation and no_slip_top were set to True.")
                Call stdout%print("        No_slip_top = True will be retained.")
                Call stdout%print("        Strict_L_Conservation will be set to False.")
                Call stdout%print("        If you wish to set Strict_L_Conservation=True, please explictly")
                Call stdout%print("        set no_slip_top to False in main_input.")      
                Call stdout%print("        ")      
            Endif          
            strict_L_conservation = .false.        
        Endif

        If (no_slip_boundaries .and. strict_L_Conservation) Then
            If (my_rank .eq. 0) Then
                Call stdout%print(" -- Warning:  Incompatible boundary conditions.")
                Call stdout%print("        Both strict_L_Conservation and no_slip_boundaries were set to True.")
                Call stdout%print("        no_slip_boundaries = True will be retained.")
                Call stdout%print("        Strict_L_Conservation will be set to False.")
                Call stdout%print("        If you wish to set Strict_L_Conservation=True, please explictly")
                Call stdout%print("        set no_slip_boundaries to False in main_input.")   
                Call stdout%print("        ")               
            Endif         
            strict_L_conservation = .false.                  
        Endif

        If (no_slip_boundaries) Then
            no_slip_top = .true.
            no_slip_bottom = .true.
        Endif
        stress_free_top    = .not. no_slip_top
        stress_free_bottom = .not. no_slip_bottom




        !/////////////////////////////////////////
        ! There are three possible magnetic boundary conditions
        ! 1.  (default)  [potential field top, potential field bottom]
        ! 2.  (dipole field bottom) [potential field top, dipole field bottom]
        ! 3.  (impose dipole field) [dipole field top, dipole field bottom]

        If (impose_dipole_field) fix_poloidalfield_top = .true.

        If (impose_dipole_field .or. dipole_field_bottom) Then
            ! Set B.C.'s on the poloidal vector
            ! potential using the value of Br at
            ! the lower boundary.

            fix_poloidalfield_bottom = .true.
            tilt_angle_radians = pi/180.0*dipole_tilt_degrees
            a = cos(tilt_angle_radians)*sqrt(4.0d0*Pi/3.0d0)
            b = sin(tilt_angle_radians)*sqrt(8.0d0*Pi/3.0d0)

            ! We use the bottom values to derive the top values
            C10_bottom = a*Br_bottom/2.0d0*(radius(N_r)**2)
            C10_top = C10_bottom*(radius(N_R)/radius(1))

            C11_bottom = b*Br_bottom/2.0d0*(radius(N_r)**2)
            C11_top = C11_bottom*(radius(N_R)/radius(1))

            C1m1_bottom = 0.0d0 ! For the time being, we don't
            C1m1_top =  0.0d0   ! worry about phasing in longitude.

        Endif

        If (adjust_dTdr_top .and. fix_dtdr_top) Then
            If (heating_type .gt. 0) Then

                ! Adjust dTdr_top so that it is consistent with internal heating and lower B.C.

                Ftop = ra_constants(10)/(four_pi*radius(1)**2) ! Flux due to internal heating

                If (fix_dtdr_bottom) Then
                    ! Account for any additional flux passing through the lower boundary
                    Fbottom=-dTdr_bottom*ref%density(N_R)*ref%temperature(N_R)*kappa(N_R)
                    Ftop = Ftop+Fbottom*(radius(N_R)/radius(1))**2
                Endif
                
                rhotk_top = ref%density(1)*ref%temperature(1)*kappa(1)
                dTdr_top = -Ftop/rhotk_top

            Endif
        Endif

        Call Generate_Boundary_Mask()

        Call Transport_Dependencies()

    End Subroutine Initialize_Boundary_Conditions

    Subroutine Finalize_Boundary_Conditions()
        Implicit None
        If (full_restart) Then
            ! Grab the boundary conditions from the checkpoint directory
            ! This will typically be used in full_restart mode.
            Call Load_BC_Mask(bc_values)  ! from checkpointing

        Endif
    End Subroutine Finalize_Boundary_Conditions

    Subroutine Generate_Boundary_Mask()
        Implicit None
        Integer :: i
        Integer :: uind, lind
        Integer :: real_ind, imag_ind
        Real*8 :: bc_val


        allocate(bc_values(2, 2, my_num_lm, n_equations))
        bc_values = zero

        uind = 1   ! Upper boundary in BC array
        lind = 2   ! Lower boundary in BC array
        real_ind = 1 ! real index in BC array
        imag_ind = 2 ! imaginary index in BC array

        !BC's are specified in physical space, but enforced in spectral space.
        !Need to multiply by sqrt(4pi) for ell=0 BCs as a result.

        If (.not. full_restart) Then

            If (fix_tvar_top .or. couple_tvar_top) Then
                if (trim(T_top_file) .eq. '__nothing__') then
                  bc_val= T_Top*sqrt(four_pi)
                  Call Set_BC(bc_val,0,0, teq,real_ind, uind)
                else  
                  Call Set_BC_from_file(T_top_file, teq, uind)
                end if
            Endif

            If (fix_tvar_bottom .or. couple_tvar_bottom) Then
                if (trim(T_bottom_file) .eq. '__nothing__') then
                  bc_val= T_bottom*sqrt(four_pi)
                  Call Set_BC(bc_val,0,0, teq,real_ind, lind)
                else
                  Call Set_BC_from_file(T_bottom_file, teq, lind)
                end if
            Endif

            If (fix_dtdr_top .or. couple_dtdr_top) Then
                if (trim(dTdr_top_file) .eq. '__nothing__') then
                  bc_val= dtdr_top*sqrt(four_pi)
                  Call Set_BC(bc_val,0,0, teq,real_ind, uind)
                else
                  Call Set_BC_from_file(dTdr_top_file, teq, uind)
                end if
            Endif

            If (fix_dtdr_bottom .or. couple_dtdr_bottom) Then
                if (trim(dTdr_bottom_file) .eq. '__nothing__') then
                  bc_val= dtdr_bottom*sqrt(four_pi)
                  Call Set_BC(bc_val,0,0, teq,real_ind, lind)
                else
                  Call Set_BC_from_file(dTdr_bottom_file, teq, lind)
                end if
            Endif

            If (fix_poloidalfield_top) Then
                if (trim(C_top_file) .eq. '__nothing__') then
                  Call Set_BC(C10_top ,1,0,ceq,real_ind,uind)
                  Call Set_BC(C11_top ,1,1,ceq,real_ind,uind)
                  Call Set_BC(C1m1_top,1,1,ceq,imag_ind,uind)
                else
                  Call Set_BC_from_file(C_top_file, ceq, uind)
                end if
            Endif        

            If (fix_poloidalfield_bottom) Then
                if (trim(C_bottom_file) .eq. '__nothing__') then
                  Call Set_BC(C10_bottom ,1,0,ceq,real_ind,lind)
                  Call Set_BC(C11_bottom ,1,1,ceq,real_ind,lind)
                  Call Set_BC(C1m1_bottom,1,1,ceq,imag_ind,lind)
                else
                  Call Set_BC_from_file(C_bottom_file, ceq, lind)
                end if
            Endif   

            do i = 1, n_active_scalars
              If (fix_chivar_a_top(i) .or. couple_chivar_a_top(i)) Then
                  if (trim(chi_a_top_file(i)) .eq. '__nothing__') then
                    bc_val= chi_a_Top(i)*sqrt(four_pi)
                    Call Set_BC(bc_val,0,0, chiaeq(i),real_ind, uind)
                  else  
                    Call Set_BC_from_file(chi_a_top_file(i), chiaeq(i), uind)
                  end if
              Endif
            end do

            do i = 1, n_passive_scalars
              If (fix_chivar_p_top(i)) Then
                  if (trim(chi_p_top_file(i)) .eq. '__nothing__') then
                    bc_val= chi_p_Top(i)*sqrt(four_pi)
                    Call Set_BC(bc_val,0,0, chipeq(i),real_ind, uind)
                  else  
                    Call Set_BC_from_file(chi_p_top_file(i), chipeq(i), uind)
                  end if
              Endif
            end do

            do i = 1, n_active_scalars
              If (fix_chivar_a_bottom(i) .or. couple_chivar_a_bottom(i)) Then
                  if (trim(chi_a_bottom_file(i)) .eq. '__nothing__') then
                    bc_val= chi_a_bottom(i)*sqrt(four_pi)
                    Call Set_BC(bc_val,0,0, chiaeq(i),real_ind, lind)
                  else
                    Call Set_BC_from_file(chi_a_bottom_file(i), chiaeq(i), lind)
                  end if
              Endif
            end do

            do i = 1, n_passive_scalars
              If (fix_chivar_p_bottom(i)) Then
                  if (trim(chi_p_bottom_file(i)) .eq. '__nothing__') then
                    bc_val= chi_p_bottom(i)*sqrt(four_pi)
                    Call Set_BC(bc_val,0,0, chipeq(i),real_ind, lind)
                  else
                    Call Set_BC_from_file(chi_p_bottom_file(i), chipeq(i), lind)
                  end if
              Endif
            end do

            do i = 1, n_active_scalars
              If (Fix_dchidr_a_Top(i) .or. couple_dchidr_a_top(i)) Then
                  if (trim(dchidr_a_top_file(i)) .eq. '__nothing__') then
                    bc_val= dchidr_a_top(i)*sqrt(four_pi)
                    Call Set_BC(bc_val,0,0, chiaeq(i),real_ind, uind)
                  else
                    Call Set_BC_from_file(dchidr_a_top_file(i), chiaeq(i), uind)
                  end if
              Endif
            end do

            do i = 1, n_passive_scalars
              If (Fix_dchidr_p_Top(i)) Then
                  if (trim(dchidr_p_top_file(i)) .eq. '__nothing__') then
                    bc_val= dchidr_p_top(i)*sqrt(four_pi)
                    Call Set_BC(bc_val,0,0, chipeq(i),real_ind, uind)
                  else
                    Call Set_BC_from_file(dchidr_p_top_file(i), chipeq(i), uind)
                  end if
              Endif
            end do

            do i = 1, n_active_scalars
              If (Fix_dchidr_a_Bottom(i) .or. couple_dchidr_a_bottom(i)) Then
                  if (trim(dchidr_a_bottom_file(i)) .eq. '__nothing__') then
                    bc_val= dchidr_a_bottom(i)*sqrt(four_pi)
                    Call Set_BC(bc_val,0,0, chiaeq(i),real_ind, lind)
                  else
                    Call Set_BC_from_file(dchidr_a_bottom_file(i), chiaeq(i), lind)
                  end if
              Endif
            end do

            do i = 1, n_passive_scalars
              If (Fix_dchidr_p_Bottom(i)) Then
                  if (trim(dchidr_p_bottom_file(i)) .eq. '__nothing__') then
                    bc_val= dchidr_p_bottom(i)*sqrt(four_pi)
                    Call Set_BC(bc_val,0,0, chipeq(i),real_ind, lind)
                  else
                    Call Set_BC_from_file(dchidr_p_bottom_file(i), chipeq(i), lind)
                  end if
              Endif
            end do

        Endif
        Call Store_BC_Mask(bc_values)  ! to checkpointing
    End Subroutine Generate_Boundary_Mask

    Subroutine Set_BC(inval,ellval,emval,eqval,imi,bound_ind)
        Implicit None
        Real*8 , Intent(In) :: inval
        Integer, Intent(In) :: ellval, emval, eqval,bound_ind, imi
        Integer :: this_l, this_m, lp, lpi
        ! Probably more efficient to have a support array / reverse lookup table
        ! In practice, this will still be short loop because
        ! There are typically only a few lm values in this config
        Do lp = 1, my_num_lm
            ! This loop is over the locally owned lm modes but the indexing arrays
            ! below are over the l and m indices of this entire radial group
            ! so we have to index into them using my_lm_min
            lpi = my_lm_min + lp - 1
            this_l = l_lm_values(lpi)
            this_m = m_lm_values(lpi)
            If ( (this_l .eq. ellval) .and. (this_m .eq. emval) ) Then
                bc_values(bound_ind , imi , lp,  eqval ) = inval
            Endif
        Enddo        
    End Subroutine Set_BC

    Subroutine Set_BC_from_file(filename,eqval,bound_ind)
        Implicit None
        Character*120, Intent(in) :: filename
        Integer, Intent(In) :: eqval,bound_ind
        integer :: fcount(3,2), r_ind
        type(SphericalBuffer) :: tempfield

        r_ind = 1
        if (bound_ind .eq. 2) then
          r_ind = n_r
        end if

        fcount(:,:) = 1
        call tempfield%init(field_count = fcount, config = 'p1b')
        call tempfield%construct('p1b')

        call read_input(filename, 1, tempfield)

        call tempfield%construct('p1a')
        call gridcp%from_spectral(tempfield%p1b, tempfield%p1a)
        call tempfield%deconstruct('p1b')

        bc_values(bound_ind,:,:,eqval) = tempfield%p1a(r_ind,:,:,1)
        call tempfield%deconstruct('p1a')

    End Subroutine Set_BC_from_file

    Subroutine Transport_Dependencies()
        Implicit None
        Real*8 :: lum_top, lum_bottom

        ! Odd reference state quantities that somehow depend on 
        ! boundary conditions can be set here.

        If (adjust_reference_heating ) Then
            !Renormalize the reference heating based on boundary fluxes
            lum_top    = 0.0d0
            lum_bottom = 0.0d0

            If ( fix_dtdr_top ) Then
                lum_top = -dtdr_top*kappa(1)*four_pi*(rmax**2)
                lum_top = lum_top*ref%density(1)*ref%temperature(1)
            Endif

            If ( fix_dtdr_bottom) Then
                lum_bottom = -dtdr_bottom*kappa(N_R)*four_pi*(rmin**2)
                lum_bottom = lum_bottom*ref%density(N_R)*ref%temperature(N_R)
            Endif

            ref%heating = ref%heating*(lum_top-lum_bottom)
            ! If the code has reached this point, c_10 was not defined in
            ! initialize_reference_heating and needs to be defined here
            ra_constants(10) = lum_top - lum_bottom
            !If something has been set inconsistenly, this will result
            ! in zero reference heating
         Endif
    End Subroutine Transport_Dependencies



    Subroutine Restore_BoundaryCondition_Defaults()
        Implicit None
        Fix_Tvar_Top    = .False.
        Fix_Tvar_Bottom = .False.
        Fix_dTdr_Top    = .False.
        Fix_dTdr_Bottom = .False.
        Fix_divrt_top   = .False.
        Fix_divt_top    = .False.
        Fix_divrfc_top  = .False.
        Fix_divfc_top   = .False.
        Fix_poloidalfield_top    = .False.
        Fix_poloidalfield_bottom = .False.
        Impose_Dipole_Field      = .False.

        T_Bottom     = 0.0d0
        T_Top        = 0.0d0
        dTdr_Top     = 0.0d0
        dTdr_Bottom  = 0.0d0
        C10_bottom = 0.0d0
        C10_top = 0.0d0
        C11_bottom = 0.0d0
        C11_top = 0.0d0
        C1m1_bottom = 0.0d0
        C1m1_top = 0.0d0
        Br_bottom = 0.0d0
        Dipole_Tilt_Degrees = 0.0d0

        T_top_file       = '__nothing__'
        T_bottom_file    = '__nothing__'
        dTdr_top_file    = '__nothing__'
        dTdr_bottom_file = '__nothing__'
        C_top_file       = '__nothing__'
        C_bottom_file    = '__nothing__'

        Strict_L_Conservation = .false.
        no_slip_boundaries = .false.
    End Subroutine Restore_BoundaryCondition_Defaults
End Module BoundaryConditions
