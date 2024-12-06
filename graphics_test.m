% graphics_test.m - Check MATLAB Graphics Setup and Renderer

clc;
disp('MATLAB Graphics and Renderer Test');
disp('----------------------------------');

try
    % Create a simple figure
    figure('Name', 'Graphics Test Figure');
    plot(1:10, '-o', 'LineWidth', 2);
    title('Graphics Test');
    xlabel('X-axis');
    ylabel('Y-axis');
    grid on;
    disp('1. Plot created successfully.');
catch ME
    disp('Error creating a plot:');
    disp(ME.message);
    return;
end

try
    % Check the renderer info
    info = rendererinfo(gca);
    disp('2. Renderer Information:');
    disp(info);
catch ME
    disp('Error retrieving renderer info:');
    disp(ME.message);
end

try
    % Check OpenGL support
    disp('3. OpenGL Information:');
    opengl_info = opengl('info');
    disp(opengl_info);
catch ME
    disp('Error retrieving OpenGL info:');
    disp(ME.message);
end

try
    % Test hardware vs. software rendering
    disp('4. Testing Renderer Modes:');
    disp('   Switching to Software OpenGL...');
    opengl software;
    figure('Name', 'Software OpenGL Test');
    plot(1:10, '-o', 'LineWidth', 2);
    title('Software OpenGL Test');
    disp('   Software OpenGL works.');

    disp('   Switching to Hardware OpenGL...');
    opengl hardware;
    figure('Name', 'Hardware OpenGL Test');
    plot(1:10, '-o', 'LineWidth', 2);
    title('Hardware OpenGL Test');
    disp('   Hardware OpenGL works.');

catch ME
    disp('Error during renderer mode tests:');
    disp(ME.message);
end

disp('----------------------------------');
disp('Graphics test completed.');
