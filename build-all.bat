@echo off

setlocal enableextensions
md "build\%COMPUTERNAME%"
endlocal


for /l %%S in (1,1,20) do (
	for /l %%E in (1,1,%%S) do (
	    set CMDLINE_DEFINES=/D "K_SYS=%%S" /D "K_ENV=%%E"
		msbuild /p:Configuration=Release
		move Release\CoffeeCode.exe build\%COMPUTERNAME%\CoffeeCode.%%S.%%E.exe
	)
)