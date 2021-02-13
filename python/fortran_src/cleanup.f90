      subroutine cleanup(did_timestep)
        ! This subroutine mimics the cleanup done in fileout_par() and
        ! fileout() (both in fileout.f) except it is called on all
        ! procs, not just proc 0.  This routine also deallocates some
        ! arrays that are allocated in runvmec.f.  This cleanup is
        ! necessary in order to read in an input file after another
        ! input file has already been read.

        ! The parameter did_timestep should be set to false when
        ! calling this routine after reading an input file but not
        ! time stepping. The parameter should be set to true when
        ! calling this routine after time stepping.
        
        ! Originally written by Matt Landreman, 2020-02-12
        
        use vmec_main
        USE vmec_params, ONLY: mscale, nscale, uminus
        use vac_persistent
        use parallel_vmec_module, ONLY: RUNVMEC_COMM_WORLD, NS_COMM,&
             FinalizeSurfaceComm, FinalizeRunVmec, grid_procs,&
             grid_size, grid_time, f3d_time, f3d_num, vgrid_time
        
        implicit none

        logical, intent(in) :: did_timestep
        integer :: istat1, ierr

        ! This barrier is here because without it, the wout file is
        ! not written for >1 proc. I don't understand why this
        ! happens, but the barrier seems to resolve it.
        ierr = 0
        call mpi_barrier(RUNVMEC_COMM_WORLD, ierr)
        
        IF (ALLOCATED(cosmu)) &
             DEALLOCATE(cosmu, sinmu, cosmum, sinmum, cosmui, cosmumi,&
             sinmui, sinmumi, cosnv, sinnv, cosnvn, sinnvn,&
             cosmui3, cosmumi3, cos01, sin01, stat=istat1)
        IF (istat1 .ne. 0) PRINT *,  "cleanup deallocate error #1"

        IF (ALLOCATED(xm)) DEALLOCATE (xm, xn, ixm, xm_nyq, xn_nyq, &
             jmin3, mscale, nscale, uminus, stat=istat1)
        IF (istat1 .ne. 0) PRINT *, "cleanup deallocate error #2"

        IF (ALLOCATED(tanu)) &
             DEALLOCATE(tanu, tanv, sinper, cosper, sinuv, cosuv, sinu,&
             cosu, sinv, cosv, sinui, cosui, cmns, csign, sinu1,&
             cosu1, sinv1, cosv1, imirr, xmpot, xnpot,&
             stat=istat1)
        IF (istat1 .ne. 0) PRINT *, "cleanup deallocate error #3"
     
        CALL free_persistent_mem 
        CALL free_mem_funct3d
        CALL free_mem_ns (.false.)
        CALL free_mem_nunv
        CALL close_all_files

        ! Now deallocate some arrays that are allocated in runvmec:
        IF(ALLOCATED(grid_procs)) THEN
           DEALLOCATE(grid_procs)
           DEALLOCATE(grid_size)
           DEALLOCATE(grid_time)
           DEALLOCATE(f3d_time)
           DEALLOCATE(f3d_num)
           IF (lfreeb) DEALLOCATE(vgrid_time)
        END IF

        if (did_timestep) then
           call FinalizeSurfaceComm(NS_COMM)
           call FinalizeRunVmec(RUNVMEC_COMM_WORLD)
        end if
        
      end subroutine cleanup
