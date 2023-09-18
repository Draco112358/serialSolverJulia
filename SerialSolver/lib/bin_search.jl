function bin_search(num, A)
    index = 0
    n = length(A)
    left = 1
    right = n
    while left <= right
        mid = ceil((left + right) / 2)
        if A[mid] == num
            index = mid
            break
        else
            if A[mid] > num
                right = mid - 1
            else
                left = mid + 1
            end
        end
    end
    return convert(Int64,index)
end