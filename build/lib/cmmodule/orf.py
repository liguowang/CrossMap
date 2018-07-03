import sys
class ORFFinder:
  """Find the longest ORF in a given sequence 
  
   "seq" is a string, if "start" is not provided any codon can be the start of 
   and ORF. If muliple ORFs have the longest length the first one encountered
   is printed 
   """
  def __init__(self, seq):
    self.seq = seq.upper()
    self.result = ("+",0,0,0,0)
    self.winner = 0
  
  def _reverse_comp(self):
    swap = {"A":"T", "T":"A", "C":"G", "G":"C","N":"N","X":"X"}
    return "".join(swap[b] for b in self.seq)[::-1]
  
  def codons(self, frame):
    """ A generator that yields DNA in one codon blocks 
    
    "frame" counts for 0. This function yelids a tuple (triplet, index) with 
    index relative to the original DNA sequence 
    """
    start = frame
    while start + 3 <= len(self.seq):
      yield (self.seq[start:start+3], start)
      start += 3 
 
  def run_one(self, frame_number, direction,start_coden, stop_coden):
    """ Search in one reading frame """
    codon_gen = self.codons(frame_number)  
    start_codens = start_coden
    stop_codens = stop_coden   
    while True:
      try: 
        c , index = codon_gen.next()
      except StopIteration:
        break 
      # Lots of conditions here: checks if we care about looking for start 
      # codon then that codon is not a stop
      if c in start_codens or not start_codens and c not in stop_codens:
        orf_start = index  # we'll return the result as 0-indexed
        end = False
        while True:
          try: 
            c, index = codon_gen.next()
          except StopIteration:
            end = True
          if c in stop_codens:
            end = True
          if end:
            orf_end = index + 3 # because index is realitve to start of codon
            L = (orf_end - orf_start)
            if L > self.winner:
              self.winner = L
              self.result = (direction, frame_number+1, orf_start, orf_end, L)
            break
    
  def longest_orf(self,direction,start_coden=['ATG'], stop_coden=['TAG','TAA','TGA']):
    if direction == "+":
      for frame in range(3):
        self.run_one(frame, direction,start_coden, stop_coden)
      return (self.result[4], self.result[1],self.seq[self.result[2]:self.result[3]]) #CDS length, coding frame, CDS sequence
    if direction == "-":
      self.seq = self._reverse_comp()
      for frame in range(3):
        self.run_one(frame, direction,start_coden, stop_coden)
      return (self.result[4], self.result[1],self.seq[self.result[2]:self.result[3]]) #CDS length, coding frame, CDS sequence

#===================
def little_test():
  seq=''
  for line in open(sys.argv[1],'r'):
    line=line.rstrip('\n\r')
    if line.startswith('>'):
  	  continue
    seq	+= line
  (l,f,s) = ORFFinder(seq).longest_orf(sys.argv[2])
  print str(l) + '\t' + str(f) + '\t' + s
  
if __name__ == "__main__":
  little_test()