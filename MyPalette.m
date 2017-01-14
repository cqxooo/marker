function theCol = MyPalette( theNumber )
while theNumber>19
    theNumber = theNumber-19;
end
theColVec = NiceColours(theNumber);

theCol = theColVec(theNumber,:);