import numpy as np


def getData(filename):
    file = open(filename, "r")
    matrix = [line.split() for line in file]
    matrix = [[float(d) for d in i] for i in matrix]
    return matrix


# zamiana macierzy PP na PD, rozłożenie układu równań na postać MACIERZxVECTOR=TARGET
def makeDualMatrix(matrix):
    eqnum = len(matrix)
    varnum = len(matrix[eqnum - 1])
    conditionnum = len(matrix) - 1
    dualmatrix = []
    vector = []
    target = []
    for i in range(conditionnum):
        target.append(matrix[i][varnum])
    for j in range(0, varnum, 1):
        dualmatrix.append([])
        for i in range(0, eqnum, 1):
            if (i == eqnum - 1):
                vector.append(matrix[i][j])
            else:
                dualmatrix[j].append(matrix[i][j])
    return dualmatrix, vector, target


# solving equations
def solveEq(eq1, eq2, sol1, sol2):
    mat = []
    vec = []
    mat.append(eq1)
    mat.append(eq2)
    vec.append(sol1)
    vec.append(sol2)
    result = np.linalg.solve(mat, vec)
    x = result[0]
    y = result[1]
    result = [x, y]
    return result


# check if the point is not yet on the list
def notHere(point, all):
    bool = True
    if (all == []):
        return bool
    else:
        for element in all:
            if (point[0] == (element[0])):
                if (point[1] == element[1]):
                    bool = False
                    return bool
            else:
                bool = True
    return bool


# check if the point is in the first quarter
def inRange(point):
    if (point[0] > 0 and point[1] > 0):
        return True


# find all points on ax X
def getAllCrossingAxXPoints(dualmatrix, vector):
    points = []
    for i in range(len(dualmatrix)):
        if (dualmatrix[i][0] != 0):
            point = [(vector[i]) / dualmatrix[i][0], 0]
        points.append(point)
    return points


# find all points on ax Y
def getAllCrossingAxYPoints(dualmatrix, vector):
    points = []
    for i in range(len(dualmatrix)):
        if (dualmatrix[i][1] != 0):
            point = [0, (vector[i]) / dualmatrix[i][1]]
        points.append(point)
    return points


# find all points of intersection of lines NOT on axes XY
def getCrossPoints(dualmatrix, vector):
    results = []
    result = []
    for i in range(len(dualmatrix)):
        for j in range(len(dualmatrix)):
            if (i != j and dualmatrix[i] != dualmatrix[j]):
                result = solveEq(dualmatrix[i], dualmatrix[j], vector[i], vector[j])
                if (notHere(result, results) and inRange(result)):
                    results.append(result)
    return (results)


# put all found points together
def setOfPoints(dualmatrix, vector):
    points = getCrossPoints(dualmatrix, vector)
    pX = getAllCrossingAxXPoints(dualmatrix, vector)
    pY = getAllCrossingAxYPoints(dualmatrix, vector)
    for point in pX:
        if (notHere(point, points)):
            points.append(point)
    for point in pY:
        if (notHere(point, points)):
            points.append(point)
    return points


# select points that meet conditions
def selectPointsInField(dualmatrix, points):
    result = []
    for i in range(len(points)):
        point = points[i]
        for j in range(len(dualmatrix)):
            value = points[i][0] * dualmatrix[j][0] + points[i][1] * dualmatrix[j][1]
            if (value < vector[j]):
                bool = False
                break
            else:
                bool = True
        if (bool == True):
            if (notHere(point, result)):
                result.append(point)
    return result


# find best point that meets conditions
def findOptimum(target, points):
    results = []
    optimumpoint = []
    for point in points:
        result = 0
        result += target[0] * point[0] + target[1] * point[1]
        results.append(result)
    minimum = min(results)
    optimumpoint = points[results.index(minimum)]
    return optimumpoint, minimum


# find lines on boundary condition
def findLinesOnBoundaryCondition(optimumPoint, dualmatrix):
    lines = []
    for i in range(len(dualmatrix)):
        if (optimumPoint[0] * dualmatrix[i][0] + optimumPoint[1] * dualmatrix[i][1] == vector[i]):
            lines.append(i)
    return lines


# znalezienie Punktu V = (x1, x2, ... , xn) realizujący optimum PP
def solvePP(matrix, lines):
    # przerabiam macierz na taką, która nie będzie miała pól które sie zerują, bo spełniają warunki silnie '<' '>'
    matrix.remove(matrix[-1])
    result = []
    target = []
    for i in range(len(matrix)):
        varnum = len(matrix[i]) - 1
        target.append(matrix[i][-1])
        primmatrix = []
        for j in range(len(matrix[i])):
            if (j in lines):
                primmatrix.append(matrix[i][j])
        result.append(primmatrix)

    A = np.array(result)
    b = np.array(target)
    x = np.linalg.solve(A, b)
    final = []
    counter = 0
    for i in range(varnum):
        if (i not in lines):
            final.append(0)
        else:
            final.append(x[counter])
            counter += 1
    return (final)


# data from file
matrix = getData("data.txt")
print(matrix)
print("\n")
# data base for PD
dualmatrix, vector, target = makeDualMatrix(matrix)
print(dualmatrix)
print("wektor wyrazów wolnych",vector)
print("Funkcja celu", target)

# all crossing points
points = setOfPoints(dualmatrix, vector)
# all points that meet the conditions
result = selectPointsInField(dualmatrix, points)
# best point in PD progrsm
optimumPointPD, maximum = findOptimum(target, result)
print("\n")
print("Punkt optymalny PD:",optimumPointPD)

lines = findLinesOnBoundaryCondition(optimumPointPD, dualmatrix)
optimumPointPP = solvePP(matrix, lines)
print("\n")
print("Lista punktów ograniczających zbiór rozwiązań dopuszczalnych dla PD:", result)
print("Punkt V realizujący optimum PP = ", optimumPointPP)
print("Wartość maksymalną: F(V)=", maximum)
