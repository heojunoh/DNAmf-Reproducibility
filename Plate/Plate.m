function [freq] = Plate(a1, a2, a3, t)
    % Create FEM model for structural modal analysis
    model = femodel(AnalysisType="structuralModal", Geometry="Plate10x10x1.stl");
    
    % Set material properties
    model.MaterialProperties = materialProperties(YoungsModulus=(1+4*a1)*1e11, ... % 1-5
                                                 PoissonsRatio=0.2+0.2*a2, ... % 0.2-0.4
                                                 MassDensity=(6+3*a3)*1e3); % 6-9
    
    % Apply boundary conditions
    model.FaceBC(1:4) = faceBC(ZDisplacement=0);
    
    % Generate mesh with max edge length
    model = generateMesh(model, Hmax=t);
    
    % Reference frequencies in Hz
    refFreqHz = [0 0 0 45.897 109.44 109.44 167.89 193.59 206.19 206.19];
    
    % Solve for natural frequencies
    maxFreq = 1.1 * refFreqHz(end) * 2 * pi;
    result = solve(model, FrequencyRange=[-0.1 maxFreq]);
    
    % Compute natural frequencies in Hz
    freqHz = result.NaturalFrequencies / (2 * pi);
    
    % Compare reference and computed frequency
    i = 4;
    tfreqHz = table(refFreqHz(i).', freqHz(i));
    tfreqHz.Properties.VariableNames = {'Reference', 'Computed'};

    freq = tfreqHz.Computed;
end
