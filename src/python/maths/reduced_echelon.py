# from https://rosettacode.org/wiki/Reduced_row_echelon_form#Python


def reduced_row_echelon_form(m):
    if not m: return
    reduced_matrix = [[m[i][j] for j in range(len(m[0]))] for i in range(len(m))]
    lead = 0
    rowCount = len(reduced_matrix)
    columnCount = len(reduced_matrix[0])
    for r in range(rowCount):
        if lead >= columnCount:
            return reduced_matrix
        i = r
        while reduced_matrix[i][lead] == 0:
            i += 1
            if i == rowCount:
                i = r
                lead += 1
                if columnCount == lead:
                    return reduced_matrix
        reduced_matrix[i], reduced_matrix[r] = reduced_matrix[r], reduced_matrix[i]
        lv = reduced_matrix[r][lead]
        reduced_matrix[r] = [mrx / float(lv) for mrx in reduced_matrix[r]]
        for i in range(rowCount):
            if i != r:
                lv = reduced_matrix[i][lead]
                reduced_matrix[i] = [iv - lv * rv for rv, iv in zip(reduced_matrix[r], reduced_matrix[i])]
        lead += 1
    return reduced_matrix

def identity(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]

def catenate(m, n):
    return [list(m[i])+list(n[i]) for i in range(len(m))]

if __name__ == '__main__':
    mtx = [
        [1, 2, -1, -4],
        [2, 3, -1, -11],
        [-2, 0, -3, 22], ]
    print(reduced_row_echelon_form(mtx))