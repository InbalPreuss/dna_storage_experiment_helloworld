def color_diff_hamming(str1, str2):
    # Calculate the Levenshtein distance between the two strings
    dist = [[0 for _ in range(len(str2) + 1)] for _ in range(len(str1) + 1)]
    for i in range(1, len(str1) + 1):
        dist[i][0] = i
    for j in range(1, len(str2) + 1):
        dist[0][j] = j
    for j in range(1, len(str2) + 1):
        for i in range(1, len(str1) + 1):
            if str1[i - 1] == str2[j - 1]:
                cost = 0
            else:
                cost = 1
            dist[i][j] = min(dist[i - 1][j] + 1, dist[i][j - 1] + 1, dist[i - 1][j - 1] + cost)

    # Color the differing letters in the original strings and return the result
    result = ""
    for i in range(len(str1)):
        if i >= len(str2) or str1[i] != str2[i]:
            result += "\033[91m" + str1[i] + "\033[0m"
        else:
            result += str1[i]
    if len(str2) > len(str1):
        result += "\033[91m" + str2[len(str1):] + "\033[0m"
    return dist[i][i], result


def color_diff_without_sub_err_considration(str1, str2):
    # Calculate the Levenshtein distance between the two strings
    m, n = len(str1), len(str2)
    dist = [[0] * (n + 1) for i in range(m + 1)]
    for i in range(1, m + 1):
        dist[i][0] = i
    for j in range(1, n + 1):
        dist[0][j] = j
    for j in range(1, n + 1):
        for i in range(1, m + 1):
            if str1[i - 1] == str2[j - 1]:
                cost = 0
            else:
                cost = 1
            dist[i][j] = min(dist[i - 1][j] + 1, dist[i][j - 1] + 1, dist[i - 1][j - 1] + cost)

    # Color the differing letters in the original strings and return the result
    result = ""
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and str1[i - 1] == str2[j - 1]:
            result = str1[i - 1] + result
            i, j = i - 1, j - 1
        else:
            if j > 0 and (i == 0 or dist[i][j - 1] <= dist[i - 1][j] and dist[i][j - 1] <= dist[i - 1][j - 1]):
                result = '\033[91m' + str2[j - 1] + '\033[0m' + result
                j = j - 1
            elif i > 0 and (j == 0 or dist[i - 1][j] <= dist[i][j - 1] and dist[i - 1][j] <= dist[i - 1][j - 1]):
                result = '\033[91m' + str1[i - 1] + '\033[0m' + result
                i = i - 1
            else:
                result = '\033[91m' + str1[i - 1] + '\033[0m' + result
                i, j = i - 1, j - 1
    return dist[m][n], result


def color_diff_levenstein(str1, str2):
    m, n = len(str1), len(str2)
    dist = [[0] * (n + 1) for i in range(m + 1)]
    for i in range(1, m + 1):
        dist[i][0] = i
    for j in range(1, n + 1):
        dist[0][j] = j
    for j in range(1, n + 1):
        for i in range(1, m + 1):
            if str1[i - 1] == str2[j - 1]:
                cost = 0
            else:
                cost = 1
            dist[i][j] = min(dist[i - 1][j] + 1, dist[i][j - 1] + 1, dist[i - 1][j - 1] + cost)

    result = ""
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and str1[i - 1] == str2[j - 1]:
            result = str1[i - 1] + result
            i, j = i - 1, j - 1
        else:
            if j > 0 and (i == 0 or dist[i][j - 1] <= dist[i - 1][j] and dist[i][j - 1] <= dist[i - 1][j - 1]):
                result = '\033[91m' + str2[j - 1] + '\033[0m' + result
                j = j - 1
            elif i > 0 and (j == 0 or dist[i - 1][j] <= dist[i][j - 1] and dist[i - 1][j] <= dist[i - 1][j - 1]):
                result = '\033[91m' + str1[i - 1] + '\033[0m' + result
                i = i - 1
            else:
                result = '\033[91m' + str2[j - 1] + '\033[0m' + result
                i, j = i - 1, j - 1

    return dist[m][n], result



if __name__ == '__main__':
    str1 = "AACCCATCACTTACCTCTCTTAACA"
    str2 = "CCACATCTACCCTTCCACTCTATCA"

    result_hamming_dist, result_hamming = color_diff_hamming(str1=str1, str2=str2)
    print(f'hamming dist: {result_hamming_dist}, hamming result: {result_hamming}')
    result_hamming_dist, result_hamming = color_diff_hamming(str1=str2, str2=str1)
    print(f'hamming dist: {result_hamming_dist}, hamming result: {result_hamming}')
    result_levenstein_dist, result_levenstein = color_diff_levenstein(str1=str1, str2=str2)
    print(f'levenstein dist {result_levenstein_dist}, levenstein result: {result_levenstein}')
    result_levenstein_dist, result_levenstein = color_diff_levenstein(str1=str2, str2=str1)
    print(f'levenstein dist {result_levenstein_dist}, levenstein result: {result_levenstein}')
