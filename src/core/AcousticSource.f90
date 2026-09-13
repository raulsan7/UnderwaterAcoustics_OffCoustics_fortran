MODULE AcousticSource

USE Kinds
USE Types
USE WindTurbine, ONLY: WindTurbine_t

IMPLICIT NONE

PRIVATE

PUBLIC :: turbine_to_source, farm_to_source


CONTAINS

SUBROUTINE turbine_to_source(turbine, config, sources)
    CLASS(WindTurbine_t)  , INTENT(IN)  :: turbine
    TYPE(Config_t)        , INTENT(IN)  :: config
    TYPE(AcousticSource_t), INTENT(OUT) :: sources


    if (.not. ALLOCATED(turbine % F) .or. .not. ALLOCATED(turbine % x)) then
        error stop "turbine_to_source(): turbine forces not computed"
    end if

    ! Copy coordinates and active frequencies
    sources % xyz  = turbine % x
    sources % Freqs = turbine % Freqs

    ! Apply force correction
    if (config % solver % method == "ImageMethod") then
        sources % force = turbine % get_impedance_corrected_force(config % env % speed_sound_base)
    else
        sources % force = turbine % F
    end if

    write(*,*) ">> [AcousticSource] Physical turbine forces translated into mathematical point dipoles"

END SUBROUTINE turbine_to_source


SUBROUTINE farm_to_source(config, sources)
    TYPE(Config_t)        , INTENT(IN)  :: config
    TYPE(AcousticSource_t), INTENT(OUT) :: sources

    return
END SUBROUTINE farm_to_source


END MODULE AcousticSource