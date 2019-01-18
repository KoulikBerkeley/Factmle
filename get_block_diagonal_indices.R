# script for deducing the vectos blk_start and blk_end from number of blocks and dimension.

library(ramify)

blk_start = linspace(1,pp - pp/blk_number + 1, blk_number);
blk_end =  linspace(pp/blk_number, pp, blk_number)
