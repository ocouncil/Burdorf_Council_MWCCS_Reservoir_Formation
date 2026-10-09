# Collapse sequences with x base differences, sort by abundance.

# Define the cutoff
diff_cutoff = 0
strip = false

# Methods
def fasta_to_hash(infile)
  f = File.open(infile, "r")
  return_hash = {}
  name = ""
  while line = f.gets do
    if line =~ /^>/
      name = line.chomp
      return_hash[name] = ""
    else
      return_hash[name] += line.chomp.upcase
    end
  end
  f.close
  return return_hash
end

def collapse_sequence_by_x_nt_difference(seq_array, cutoff)
  new_seq_freq = {}
  seq_freq = count(seq_array)
  if seq_freq.size == 1
    new_seq_freq = seq_freq
  else
    uniq_seq = seq_freq.keys
    unique_seq_pair = uniq_seq.combination(2)
    dupli_seq = []
    unique_seq_pair.each do |pair|
      seq1, seq2 = pair
      diff = compare_two_seq(seq1, seq2)
      if diff <= cutoff
        freq1 = seq_freq[seq1]
        freq2 = seq_freq[seq2]
        freq1 >= freq2 ? dupli_seq << seq2 : dupli_seq << seq1
      end
    end
    seq_freq.each do |seq, freq|
      new_seq_freq[seq] = freq unless dupli_seq.include?(seq)
    end
  end
  return new_seq_freq
end

def compare_two_seq(seq1 = "", seq2 = "")
  length = seq1.size
  diff = 0
  (0..(length - 1)).each do |position|
    nt1 = seq1[position]
    nt2 = seq2[position]
    diff += 1 unless nt1 == nt2
  end
  return diff
end

def count(array)
  hash = Hash.new(0)
  array.each { |element| hash[element] += 1 }
  return hash
end

# Main
indir = ARGV[0]

# Create centralized output directories
collapse_dir = "centralized_collapse_output"
aln_dir = "centralized_aln_output"
summary_dir = "centralized_summary_output"

[collapse_dir, aln_dir, summary_dir].each do |dir|
  Dir.mkdir(dir) unless File.directory?(dir)
end

# files in subdirectories
f = Dir.glob("#{indir}/**/*").select { |file| File.file?(file) }

f.each do |infile|
  base = File.basename(infile, ".*")
  ext = File.extname(infile)

  align_file = "#{aln_dir}/#{base}.fa"
  summary_file = "#{summary_dir}/#{base}.csv"
  outfile = "#{collapse_dir}/#{base}_collapse#{diff_cutoff}#{ext}"

  # Run alignment
  print `muscle -in #{infile} -out #{align_file} -quiet -maxiters 2`

  out = File.open(outfile, "w")
  summary_out = File.open(summary_file, "w")
  summary_out.puts "sequence#,count"

  seqs = fasta_to_hash(align_file)
  seqs = gap_strip(seqs) if strip

  filter_seq_freq = collapse_sequence_by_x_nt_difference(seqs.values, diff_cutoff)
  n = 1
  filter_seq_freq.each do |k, v|
    name = ">#{base}_#{n}_#{v}"
    out.puts name
    out.puts k.tr("-", "")
    summary_out.puts "#{n},#{v}"
    n += 1
  end

  out.close
  summary_out.close
end
