#maksAPOBECmutations.rb
#Given a reference sequence, look for G-to-A mutations in the DNA sequences. When an “A” occurs in a position where a “G” is in the reference sequence, and the subsequent base is either a G or an A, replace the “A” mutation with “N.”

#to run in the terminal command line type:
#ruby maskAPOBECmutations.rb  nameOfFileToBeMasked.fas


inputFile = ARGV            #ARGV is a list of your input files, in this case nameOfFileToBeMasked.fas
source = inputFile[0]       #inputFile[0] picks the first (and in this case only) file that you used in the command  line


outFile = source + "_MASKED" + ".fas"  #this is the name of the masked file

writeOut = File.open(outFile, 'w')  #opens the masked file, 'w' means you can write to this file


infile = File.open(source, 'r') #opens the unmasked file, 'r' means you can read it but not change it

line_count = 1  #count through each line of the unmasked file
reference = ""  #this initializes the reference sequence, but without any info


while line = infile.gets and line.strip != nil  #read in one line of the file at a time
    if line_count == 1  #if it's the first line
        writeOut.write(line) #then write it unchanged to the output file
        elsif line_count  == 2 #if it's the second line
        reference = line  #then it's the reference sequence
        refLength = line.length()  #how long is the reference sequence?
        writeOut.write(line)  # write the unchanged reference to the output file
        elsif line_count % 2 == 1 # if it's an odd numbered line after the first line
        writeOut.write(line) #it's the name of the sequence, write it to the output file
        
        elsif line_count % 2 == 0  #if it's an even numbered line after the second line
        maskedLine = ""  #start a new sequence that will be the masked sequence
        for idx in 0..(refLength-1) # idx will count through each base of the reference and sequence
            if line.length != refLength #if the sequence isn't the same length as the reference sequence
                puts "mismatched lengths" #sound an alarm
            end
            lineChr = line[idx]  #look at the base at position idx
            if reference[idx] == "G" and lineChr == "A" and (reference[idx + 1] == "A" or reference[idx + 1] == "G")  #if it's an "A" and the reference is a "G" and the next letter is an "A" or "G"
                maskedLine += "N" #add an "N" to the masked line
                else #otherwise
                maskedLine += lineChr #just add the base to the masked line
            end
        end
        writeOut.write(maskedLine) #after going through each base in the line, write the masked line to the output file
        else
        puts line_count #this line of code should never run.  If the line_count gets printed to your terminal, you know something went wrong.
        
    end
    line_count += 1 #increase the line_count by one, and read through the while loop until there are no lines left
end


