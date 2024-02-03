function ptb3 = isptb3()
	version = PsychtoolboxVersion;
	ptb3 = ischar(version) && version(1) == '3';
end