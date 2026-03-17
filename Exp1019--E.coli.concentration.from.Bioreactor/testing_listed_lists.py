tip_col_cycles = [1, 2, 3] # how many col of tips will this require? 3.
tip_row_cycles = [["A","B"], ["C","D"], ["E","F"], ["G","H"]] # how do we group the tips? A+B, next to each other
channel_cycles = [[0,4], [1,5], [2,6], [3,7]]
ROW_PAIRS = [["A","B","C","D"], ["E","F","G","H"]]
for col in range(1,13):
    cycle_idx = (col - 1) // 4          # 0 for cols 1-4, 1 for 5-8, 2 for 9-12
    channel_idx = (col -1) % 4      # which set of channels to use
    print (channel_cycles[channel_idx])
    row_idx   = (col - 1) % 4   
    tipcol    = tip_col_cycles[cycle_idx]
    r1, r2    = tip_row_cycles[row_idx]
    tip_slice = f"{r1}{tipcol}:{r2}{tipcol}"  # e.g., "A1
