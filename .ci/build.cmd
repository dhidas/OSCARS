@echo off
:: To build extensions for 64 bit Python 3, we need to configure environment
:: variables to use the MSVC 2010 C++ compilers from GRMSDKX_EN_DVD.iso of:
:: MS Windows SDK for Windows 7 and .NET Framework 4
::
:: More details at:
:: https://github.com/cython/cython/wiki/CythonExtensionsOnWindows

IF "%DISTUTILS_USE_SDK%"=="1" (
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
) ELSE (
    ECHO Using default MSVC build environment
)

CALL %*
