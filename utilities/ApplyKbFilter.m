function ApplyKbFilter()
	WaitSecs(0.5);
	FlushEvents('keyDown');
	if isptb3
		[keyIsDown, when, keyCode] = KbCheck;
		DisableKeysForKbCheck(find(keyCode));
	end
end
