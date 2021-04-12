# The technical detail of MALBAC-DT

## Primer A / B
The experiment process:
1. Add primer A.
2. Reverse transcription (RT). (After that, the excessive primer A need to be removed.)
3. Add primer B + digestive enzyme. (In low temperature)
4. Digest excessive primer A and primer B. (The product of RT containing primer A will not be digested.)
5. PCR.

The content of primer B can reflect on the content of excessive primer A:  
* If digest efficiency is extremely high:  
Little primer A or B will be left.  
Most sequencing reads are from cDNA.  
The ratio A/(A+B) will be nearly 100%. (Only the Primer A linked to cDNA were sequenced.)

* If digest efficiency is extremely low:  
Many primer A or B will be left.  
The sequencing reads will be overwhelmed by primer A and B.  
The ratio A/(A+B) will be nearly 50%.
