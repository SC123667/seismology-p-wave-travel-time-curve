clear;
clc;
close all;

% =========================
% P 波走时曲线作业
% 说明：
% 1. 读取同文件夹里的 prem.txt
% 2. 用 PREM 的深度和 P 波速度做地球扁平化变换
% 3. 用 201 个射线参数做简化射线追踪
% 4. 画出题目要求的两张图
% =========================

earthRadius = 6371;                    % 地球半径，单位 km
rayParameter = linspace(0.0017, 0.1128, 201)';   % 题目给定的 201 个射线参数
flatStepKm = 5;                        % 平地球模型里的采样步长，单位 km
reductionVelocity = 0.1;               % 约化速度，单位 degree/s

thisFile = mfilename('fullpath');
if isempty(thisFile)
    scriptFolder = pwd;
else
    scriptFolder = fileparts(thisFile);
end

premFile = fullfile(scriptFolder, 'prem.txt');

if ~isfile(premFile)
    error('没有找到 prem.txt，请把数据文件和本脚本放在同一个文件夹中。');
end

try
    M = readmatrix(premFile, 'FileType', 'text');
catch
    M = dlmread(premFile);
end

if size(M, 2) < 3
    error('prem.txt 至少需要 3 列数据。第 1 列是深度，第 3 列是 P 波速度。');
end

depth0 = M(:, 1);
vp0 = M(:, 3);

[depth0, sortIndex] = sort(depth0, 'ascend');
vp0 = vp0(sortIndex);

% 课件里提到地心位置做扁平化会发散，所以最后一个点不要正好取到地心
if depth0(end) >= earthRadius
    depth0(end) = earthRadius - 1;
end

% 1. 地球扁平化变换
radius0 = earthRadius - depth0;
zFlat0 = earthRadius .* log(earthRadius ./ radius0);
vpFlat0 = (earthRadius ./ radius0) .* vp0;

% 2. 为了积分更稳定，把平地球模型再加密一点
[zFlat, vpFlat, depthFine] = buildFineFlatModel(zFlat0, vpFlat0, earthRadius, flatStepKm);

numRay = length(rayParameter);
distanceDegree = nan(numRay, 1);
travelTimeSecond = nan(numRay, 1);
reducedTimeSecond = nan(numRay, 1);
turnDepth = nan(numRay, 1);
phaseName = cell(numRay, 1);
branchName = cell(numRay, 1);

% 3. 对每个射线参数计算走时和距离
for i = 1:numRay
    p = rayParameter(i);

    [xOneWay, tOneWay, dTurn, branchText, isOk] = traceOneRay(p, zFlat, vpFlat, depthFine);

    if isOk
        distanceDegree(i) = 2 * xOneWay / earthRadius * 180 / pi;
        travelTimeSecond(i) = 2 * tOneWay;
        reducedTimeSecond(i) = travelTimeSecond(i) - distanceDegree(i) / reductionVelocity;
        turnDepth(i) = dTurn;
        phaseName{i} = classifyPhase(dTurn);
        branchName{i} = branchText;
    else
        phaseName{i} = 'no result';
        branchName{i} = 'no result';
    end
end

% 4. 保存结果表
resultTable = table(rayParameter, distanceDegree, travelTimeSecond, reducedTimeSecond, ...
    turnDepth, phaseName, branchName, ...
    'VariableNames', {'p_s_per_km', 'X_degree', 'T_second', 'T_reduced_second', ...
    'turn_depth_km', 'phase_name', 'branch_name'});

writetable(resultTable, fullfile(scriptFolder, 'P波走时结果表.csv'));

% 5. 图 1：X = 0 到 180°，T = 0 到 25 min
fig1 = figure('Color', 'w', 'Position', [120, 100, 980, 620]);
hold on;
grid on;
box on;

phaseList = {'P', 'PKP', 'PKIKP'};
phaseColor = [0.10 0.45 0.85;
              0.85 0.33 0.10;
              0.20 0.65 0.20];

for k = 1:length(phaseList)
    nowPhase = phaseList{k};
    mask = strcmp(phaseName, nowPhase) & ~isnan(distanceDegree) & ~isnan(travelTimeSecond);
    xData = distanceDegree(mask);
    yData = travelTimeSecond(mask) / 60;

    [xLine, yLine] = makeBrokenLine(xData, yData, 6, 2.5);
    plot(xLine, yLine, '-', 'LineWidth', 1.6, 'Color', phaseColor(k, :), 'DisplayName', nowPhase);
    plot(xData, yData, '.', 'MarkerSize', 11, 'Color', phaseColor(k, :), 'HandleVisibility', 'off');
end

xlim([0 180]);
ylim([0 25]);
xlabel('震中距 X (degree)');
ylabel('走时 T (minute)');
title('PREM 模型下的 P 波走时曲线');
legend('Location', 'northwest');

saveFigureSimple(fig1, fullfile(scriptFolder, '图1_P波走时曲线_0到180度.png'));

% 6. 图 2：X = 10 到 35°，T = 50 到 100 s，约化速度 0.1 degree/s
fig2 = figure('Color', 'w', 'Position', [160, 120, 980, 620]);
hold on;
grid on;
box on;

maskWindow = distanceDegree >= 10 & distanceDegree <= 35 & ~isnan(reducedTimeSecond);

mask1 = maskWindow & turnDepth < 400;
mask2 = maskWindow & abs(turnDepth - 400) <= 0.5;
mask3 = maskWindow & turnDepth > 400 & turnDepth < 670;
mask4 = maskWindow & abs(turnDepth - 670) <= 0.5;
mask5 = maskWindow & turnDepth > 670;

plotBranch(mask1, distanceDegree, reducedTimeSecond, [0.00 0.45 0.74], '转折深度 < 400 km');
plotBranch(mask2, distanceDegree, reducedTimeSecond, [0.85 0.33 0.10], '在 400 km 转折');
plotBranch(mask3, distanceDegree, reducedTimeSecond, [0.47 0.67 0.19], '转折深度在 400-670 km');
plotBranch(mask4, distanceDegree, reducedTimeSecond, [0.49 0.18 0.56], '在 670 km 转折');
plotBranch(mask5, distanceDegree, reducedTimeSecond, [0.93 0.69 0.13], '转折深度 > 670 km');

xlim([10 35]);
ylim([50 100]);
xlabel('震中距 X (degree)');
ylabel('约化走时 T_r (second)');
title('上地幔不连续面导致的走时曲线三重分支');
legend('Location', 'best');

saveFigureSimple(fig2, fullfile(scriptFolder, '图2_P波走时曲线_10到35度_三重分支.png'));

disp('计算完成。');
disp('已经生成：');
disp('1. 图1_P波走时曲线_0到180度.png');
disp('2. 图2_P波走时曲线_10到35度_三重分支.png');
disp('3. P波走时结果表.csv');


function [zFine, vFine, depthFine] = buildFineFlatModel(zNode, vNode, earthRadius, stepFlatKm)
zFine = [];
vFine = [];
tol = 1e-10;

for i = 1:length(zNode) - 1
    z1 = zNode(i);
    z2 = zNode(i + 1);
    v1 = vNode(i);
    v2 = vNode(i + 1);

    if abs(z2 - z1) < tol
        if isempty(zFine)
            zFine = [zFine; z1];
            vFine = [vFine; v1];
        elseif abs(zFine(end) - z1) > tol || abs(vFine(end) - v1) > tol
            zFine = [zFine; z1];
            vFine = [vFine; v1];
        end

        if abs(zFine(end) - z2) > tol || abs(vFine(end) - v2) > tol
            zFine = [zFine; z2];
            vFine = [vFine; v2];
        end
    else
        pointNumber = max(2, ceil((z2 - z1) / stepFlatKm) + 1);
        localZ = linspace(z1, z2, pointNumber)';
        localV = linspace(v1, v2, pointNumber)';

        if ~isempty(zFine)
            if abs(zFine(end) - localZ(1)) < tol && abs(vFine(end) - localV(1)) < tol
                localZ(1) = [];
                localV(1) = [];
            end
        end

        zFine = [zFine; localZ];
        vFine = [vFine; localV];
    end
end

radiusFine = earthRadius .* exp(-zFine ./ earthRadius);
depthFine = earthRadius - radiusFine;
end


function [xOneWay, tOneWay, dTurn, branchText, isOk] = traceOneRay(p, zFlat, vpFlat, depthFine)
xOneWay = 0;
tOneWay = 0;
dTurn = nan;
branchText = '';
isOk = false;

tol = 1e-10;

for i = 1:length(zFlat) - 1
    z1 = zFlat(i);
    z2 = zFlat(i + 1);
    v1 = vpFlat(i);
    v2 = vpFlat(i + 1);
    d1 = depthFine(i);

    % 先判断是不是不连续面
    if abs(z2 - z1) < tol
        s1 = 1 / v1;
        s2 = 1 / v2;

        if s2 < p && p <= s1 + tol
            dTurn = d1;
            branchText = classifyUpperMantleBranch(dTurn);
            isOk = true;
            return;
        else
            continue;
        end
    end

    f1 = p * v1 - 1;
    f2 = p * v2 - 1;

    if f2 <= tol
        [dx, dt] = integratePiece(z1, z2, v1, v2, p, 4);
        xOneWay = xOneWay + dx;
        tOneWay = tOneWay + dt;
    elseif f1 <= tol && f2 > tol
        if abs(v2 - v1) < tol
            frac = 0;
        else
            frac = (1 / p - v1) / (v2 - v1);
        end

        frac = max(0, min(1, frac));

        zTurn = z1 + frac * (z2 - z1);
        vTurn = v1 + frac * (v2 - v1);
        dTurn = depthFine(i) + frac * (depthFine(i + 1) - depthFine(i));

        [dx, dt] = integratePiece(z1, zTurn, v1, vTurn, p, 80);
        xOneWay = xOneWay + dx;
        tOneWay = tOneWay + dt;

        branchText = classifyUpperMantleBranch(dTurn);
        isOk = true;
        return;
    else
        dTurn = d1;
        branchText = classifyUpperMantleBranch(dTurn);
        isOk = true;
        return;
    end
end
end


function [dx, dt] = integratePiece(z1, z2, v1, v2, p, partNumber)
dx = 0;
dt = 0;

if z2 <= z1
    return;
end

localZ = linspace(z1, z2, partNumber + 1);
localV = linspace(v1, v2, partNumber + 1);

for k = 1:partNumber
    dz = localZ(k + 1) - localZ(k);
    vmid = (localV(k) + localV(k + 1)) / 2;

    temp = 1 - (p * vmid)^2;
    if temp < 1e-12
        temp = 1e-12;
    end

    dx = dx + dz * (p * vmid) / sqrt(temp);
    dt = dt + dz / (vmid * sqrt(temp));
end
end


function phaseText = classifyPhase(dTurn)
if isnan(dTurn)
    phaseText = 'unknown';
elseif dTurn < 2891
    phaseText = 'P';
elseif dTurn < 5150
    phaseText = 'PKP';
else
    phaseText = 'PKIKP';
end
end


function branchText = classifyUpperMantleBranch(dTurn)
if isnan(dTurn)
    branchText = 'unknown';
elseif dTurn < 400
    branchText = 'turn above 400 km';
elseif abs(dTurn - 400) <= 0.5
    branchText = 'turn at 400 km';
elseif dTurn > 400 && dTurn < 670
    branchText = 'turn between 400 and 670 km';
elseif abs(dTurn - 670) <= 0.5
    branchText = 'turn at 670 km';
else
    branchText = 'turn below 670 km';
end
end


function [xLine, yLine] = makeBrokenLine(xData, yData, xGapLimit, yGapLimit)
xLine = [];
yLine = [];

if isempty(xData)
    return;
end

for i = 1:length(xData)
    xLine = [xLine; xData(i)];
    yLine = [yLine; yData(i)];

    if i < length(xData)
        if abs(xData(i + 1) - xData(i)) > xGapLimit || abs(yData(i + 1) - yData(i)) > yGapLimit
            xLine = [xLine; nan];
            yLine = [yLine; nan];
        end
    end
end
end


function plotBranch(mask, xAll, yAll, colorValue, labelText)
xData = xAll(mask);
yData = yAll(mask);

if isempty(xData)
    plot(nan, nan, '-', 'LineWidth', 1.6, 'Color', colorValue, 'DisplayName', labelText);
    return;
end

[xLine, yLine] = makeBrokenLine(xData, yData, 1.5, 3);
plot(xLine, yLine, '-', 'LineWidth', 1.6, 'Color', colorValue, 'DisplayName', labelText);
plot(xData, yData, '.', 'MarkerSize', 12, 'Color', colorValue, 'HandleVisibility', 'off');
end


function saveFigureSimple(figHandle, fileName)
try
    exportgraphics(figHandle, fileName, 'Resolution', 300);
catch
    saveas(figHandle, fileName);
end
end
