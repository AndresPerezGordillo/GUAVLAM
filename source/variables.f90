module variables

    implicit none
    public

    real(8), dimension(:,:), allocatable :: nodes,          &
                                            control_points, &
                                            vel_cp,         &
                                            normals,        &
                                            A, A_copy,      &
                                            G,              &
                                            gam,            &
                                            plate_m_vel,    &
                                            plate_vel_jump, &
                                            plate_cp,       &
                                            wake_gam,       &
                                            wake_G,         &
                                            wake_V,         &
                                            wake_XYZ

    real(8), dimension(3) :: v_inf

    real(8) :: alfa, deltaT, span_dim, chord_dim, AR, angle_step
    real(8) :: velocity, velocity_dim
    real(8) :: L_c, V_c, span, chord

    integer, dimension(:,:), allocatable :: NB,             &
                                            panels,         &
                                            wake_panels,    &
                                            wake_nb

    integer, dimension(:), allocatable :: te_index,         &
                                          te_panels,        &
                                          w_index

    integer :: nn, np, nwn, nwp, ds, dc, disc_level, angle_sweep, i, numTSteps, t, an

    character(len=100) :: numchar

    logical :: debug_flag

end module variables
