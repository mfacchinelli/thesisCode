function returnValue = performRecursiveInterpolationStep( ...
    independentValues, ...
    dependentData, ...
    currentVariable, ...
    independentValuesToInterpolate, ...
    currentArrayIndices, ...
    nearestLowerIndices, ...
    NumberOfDimensions )

%...Calculate fractions of data points above and below independent
%...variable value to be added to interpolated value.
upperFraction = ( independentValuesToInterpolate( currentVariable ) - ...
    independentValues{ currentVariable }( nearestLowerIndices( currentVariable ) ) ) / ...
    ( independentValues{ currentVariable }( nearestLowerIndices( currentVariable ) + 1 ) - ...
    independentValues{ currentVariable }( nearestLowerIndices( currentVariable ) ) );
lowerFraction = -( independentValuesToInterpolate( currentVariable ) - ...
    independentValues{ currentVariable }( nearestLowerIndices( currentVariable ) + 1 ) ) / ...
    ( independentValues{ currentVariable }( nearestLowerIndices( currentVariable ) + 1 ) - ...
    independentValues{ currentVariable }( nearestLowerIndices( currentVariable ) ) );

%...If at top dimension, call dependent variable data.
if ( currentVariable == NumberOfDimensions )
    currentArrayIndices( NumberOfDimensions ) = nearestLowerIndices( currentVariable );
    lowerContribution = dependentData( currentArrayIndices( 1 ), currentArrayIndices( 2 ) );
    currentArrayIndices( NumberOfDimensions ) = nearestLowerIndices( currentVariable ) + 1;
    upperContribution = dependentData( currentArrayIndices( 1 ), currentArrayIndices( 2 ) );
    
    %...If at lower dimension, update currentArrayIndices and call function with
    %...currentVariable++.
else
    currentArrayIndices( currentVariable ) = nearestLowerIndices( currentVariable );
    lowerContribution = performRecursiveInterpolationStep( ...
        independentValues, dependentData, ...
        currentVariable + 1, independentValuesToInterpolate, ...
        currentArrayIndices, nearestLowerIndices, NumberOfDimensions );
    currentArrayIndices( currentVariable ) = nearestLowerIndices( currentVariable ) + 1;
    upperContribution = performRecursiveInterpolationStep( ...
        independentValues, dependentData, ...
        currentVariable + 1, independentValuesToInterpolate, ...
        currentArrayIndices, nearestLowerIndices, NumberOfDimensions );
end

%...Return interpolated value.
returnValue = upperFraction * upperContribution + ...
    lowerFraction * lowerContribution;
end