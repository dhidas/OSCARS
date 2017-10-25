@ECHO OFF

SET COMMAND_TO_RUN=%*
SET WIN_SDK_ROOT=C:\Program Files\Microsoft SDKs\Windows
SET WIN_WDK=c:\Program Files (x86)\Windows Kits\10\Include\wdf

:: Based on the Python version, determine what SDK version to use, and whether
:: to set the SDK for 64-bit.
ECHO hi1
ECHO MPV %MAJOR_PYTHON_VERSION%
IF %MAJOR_PYTHON_VERSION% == 2 (
    ECHO py2
    SET WINDOWS_SDK_VERSION="v7.0"
    SET SET_SDK_64=Y
) ELSE (
    ECHO elsepy3orgrate
    IF %MAJOR_PYTHON_VERSION% == 3 (
        ECHO py3
        SET WINDOWS_SDK_VERSION="v7.1"
        IF %MINOR_PYTHON_VERSION% LEQ 4 (
            ECHO minorlessthan4
            SET SET_SDK_64=Y
        ) ELSE (
            ECHO minorNOTlessthan4
            SET SET_SDK_64=N
            IF EXIST "%WIN_WDK%" (
                ECHO WIN_WDK exists
                :: See: https://connect.microsoft.com/VisualStudio/feedback/details/1610302/
                REN "%WIN_WDK%" 0wdf
            )
        )
    ) ELSE (
        ECHO Unsupported Python version: "%MAJOR_PYTHON_VERSION%"
        EXIT 1
    )
)

ECHO pasthefirstset
IF %ARCH% == 64 (
    IF %SET_SDK_64% == Y (
        ECHO Configuring Windows SDK %WINDOWS_SDK_VERSION% for Python %MAJOR_PYTHON_VERSION% on a 64 bit architecture
        SET DISTUTILS_USE_SDK=1
        SET MSSdk=1
        "%WIN_SDK_ROOT%\%WINDOWS_SDK_VERSION%\Setup\WindowsSdkVer.exe" -q -version:%WINDOWS_SDK_VERSION%
        "%WIN_SDK_ROOT%\%WINDOWS_SDK_VERSION%\Bin\SetEnv.cmd" /x64 /release
        ECHO Executing: %COMMAND_TO_RUN%
        call %COMMAND_TO_RUN% || EXIT 1
    ) ELSE (
        ECHO Using default MSVC build environment for 64 bit architecture
        ECHO Executing: %COMMAND_TO_RUN%
        call %COMMAND_TO_RUN% || EXIT 1
    )
) ELSE (
    ECHO Using default MSVC build environment for 32 bit architecture
    ECHO Executing: %COMMAND_TO_RUN%
    call %COMMAND_TO_RUN% || EXIT 1
)
