test_that("Input arguments are handled correctly.", {
  
  ##################################
  # test 1
  # make sure triple 0s are skipped
  t1_seq = strsplit("ATGATGATG", "")[[1]]
  t1_path =strsplit("111000112111", "")[[1]]
  t1_expected_seq = "ATGATATG"
  
  t1_frame_dat = list(framed_seq = "ATGATGATG",
                      framed_vec = t1_seq,
                      trimmed_seq = t1_seq,
                      removed_lead = NULL,
                      removed_end = NULL,
                      seq_start = 1,
                      seq_end = length(t1_seq),
                      path_start = 1 ,
                      path_end = length(t1_path),
                      front = NULL)
  
  t1_adj_out = adj_seq(t1_frame_dat, t1_path, censor_length = 0)
  t1_adj_seq = paste(t1_adj_out[[1]], collapse = "")
  #t1_adj_seq == t1_expected_seq
  expect_equal(t1_adj_seq ,  t1_expected_seq)

  #test where the censor extends past the end of the sequence
  t1_expected_seq2 = "A-------"
  
  t1_adj_out = adj_seq(t1_frame_dat, t1_path, censor_length = 4)
  t1_adj_seq = paste(t1_adj_out[[1]], collapse = "")
  #t1_adj_seq == t1_expected_seq2
  expect_equal(t1_adj_seq ,  t1_expected_seq2)
  
  
  
  ##################################
  #test 2
  #make sure inserts are handled in the correct locaiton
  t2_seq = strsplit("AAAATGATGATG", "")[[1]]
  t2_path =strsplit("111101111110", "")[[1]]
  t2_expected_seq = "AAAAT-GATGATG-"
  
  t2_frame_dat = list(framed_seq = "AAAATGATGATG",
                      framed_vec = t2_seq,
                      trimmed_seq = t2_seq,
                      removed_lead = NULL,
                      removed_end = NULL,
                      seq_start = 1,
                      seq_end = length(t2_seq),
                      path_start = 1 ,
                      path_end = length(t2_path),
                      front = NULL)
  
  t2_adj_out = adj_seq(t2_frame_dat, t2_path, censor_length = 0)
  t2_adj_seq = paste(t2_adj_out[[1]], collapse = "")
  #t2_adj_seq == t2_expected_seq
  expect_equal(t2_adj_seq ,  t2_expected_seq)
  
  t2_expected_seq2 = "AAA-----TGA---"
  # org:                   ^       ^
  #flanks:               ^^ ^^   ^^
  t2_adj_out = adj_seq(t2_frame_dat, t2_path, censor_length = 2)
  t2_adj_seq = paste(t2_adj_out[[1]], collapse = "")
  #t2_adj_seq == t2_expected_seq2
  expect_equal(t2_adj_seq ,  t2_expected_seq)
  
  
  ##################################
  #test 3
  # make sure deletes are handled in the correct location
  t3_seq = strsplit("ATGATGATG", "")[[1]]
  t3_path =strsplit("111211111", "")[[1]]
  t3_expected_seq = "ATGTGATG"
  
  t3_frame_dat = list(framed_seq = "ATGATGATG",
                      framed_vec = t3_seq,
                      trimmed_seq = t3_seq,
                      removed_lead = NULL,
                      removed_end = NULL,
                      seq_start = 1,
                      seq_end = length(t3_seq),
                      path_start = 1 ,
                      path_end = length(t3_path),
                      front = NULL)
  
  t3_adj_out = adj_seq(t3_frame_dat, t3_path, censor_length = 0)
  t3_adj_seq = paste(t3_adj_out[[1]], collapse = "")
  #t3_adj_seq == t3_expected_seq
  expect_equal(t3_adj_seq ,  t3_expected_seq)
  
  t3_expected_seq2 = "AT--GATG"
  
  t3_adj_out = adj_seq(t3_frame_dat, t3_path, censor_length = 1)
  t3_adj_seq = paste(t3_adj_out[[1]], collapse = "")
  #t3_adj_seq == t3_expected_seq2
  expect_equal(t3_adj_seq ,  t3_expected_seq2)
  
  })