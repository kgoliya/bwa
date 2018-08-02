Notes:
    pwrite() is best used with the full buffer. Make sure to wait enough time for the pwrite to finish. AWS provides fsync for this purpose but
    it is not faitfully representing completion of writes to the fpga memroy (either to dram or to on-chip memory)
