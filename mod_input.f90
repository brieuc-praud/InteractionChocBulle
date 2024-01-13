Module mod_input

    Use mod_parameters

    Implicit None

Contains
    Subroutine read_parameters(parameters)
        ! --- InOut
        Character(len=*), Intent(In) :: parameters
        ! --- Locals
        Character(len=100) :: buffer
        Character(len=50) :: label
        Integer :: pos
        Integer :: ios = 0
        Integer :: line_number = 0

        Open(FILEINPUT, File=TRIM(ADJUSTL(parameters)))

        Do While (ios == 0)
            Read(FILEINPUT, '(A)', IOstat=ios) buffer
            If (ios == 0) Then
                line_number = line_number + 1

                pos = SCAN(buffer, ' ')
                label = buffer(1:pos)
                buffer = buffer(pos+1:)

                Select Case (label)
                Case ('number_cells_x')
                    Read(buffer, *, IOstat=ios) imax
                Case ('number_cells_y')
                    Read(buffer, *, IOstat=ios) jmax
                Case ('time_max')
                    Read(buffer, *, IOstat=ios) time_max
                Case ('error_norm')
                    Read(buffer, *, IOstat=ios) norm_str
                Case ('output_modulo')
                    Read(buffer, *, IOstat=ios) output_modulo
                Case ('case')
                    Read(buffer, *, IOstat=ios) case_name
                Case ('adiabatic_index')
                    Read(buffer, *, IOstat=ios) gammagp
                Case ('CFL')
                    Read(buffer, *, IOstat=ios) cfl
                Case ('time_scheme_order')
                    Read(buffer, *, IOstat=ios) num_scheme%time_scheme_order
                Case ('flux')
                    Read(buffer, *, IOstat=ios) num_scheme%space_scheme_name
                Case ('space_scheme_order') 
                    Read(buffer, *, IOstat=ios) num_scheme%space_scheme_order
                Case ('MUSCL_limiter')
                    Read(buffer, *, IOstat=ios) num_scheme%MUSCL_limiter
                Case ('generalised_minmod_parameter')
                    Read(buffer, *, IOstat=ios) num_scheme%MUSCL_generalised_minmod_parameter
                Case ('','#')
                    ! Do nothing if it is an empty line or a comment
                Case Default
                    If ( label(1:1) /= '#' ) Then ! Special case where there is no space after '#'
                        Write(STDOUT,*) "Invalid label ", label, " at line ", line_number, " (skipping)"
                    End If
                End Select
            End If
        End Do
        Close(FILEINPUT)
    End Subroutine read_parameters


End Module mod_input
