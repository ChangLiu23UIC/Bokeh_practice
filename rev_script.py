import re
import random

"""Define a function to reverse the sequence bsed on different mode chosen
    1 is whole sequence reverse
     2 is reverse the sequence at trypsin cuts
      3 is random shuffle of the sequence
      The output will have "rev_" as prefix and in a fasta format"""

def rev_string(input: str, mode: int):
    # three modes of
    match mode:
        # reverse the whole string
        case 1:
            # generate the information from the string
            name = input.split("\n")[0] + "\n"
            # reverse the sequence
            rev_sequence = "".join(input.split("\n")[1:])[::-1]
            # format it into the result
            output = ">rev_" + name + rev_sequence + "\n"
            return output
        # reverse the sequence at trypsin cuts
        case 2:
            name = input.split("\n")[0] + "\n"
            sequence = "".join(input.split("\n")[1:])
            # Use regular expression to find the cutting site
            cut_end_pos_list = ([i.start()
                                for i in re.finditer(r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))', sequence)]
                                + [(len(sequence) - 1)])
            cut_start_pos_list = [0] + [i + 1 for i in cut_end_pos_list]
            cut_start_pos_list.pop()
            # Join the cut sequence by reversing each fragment
            rev_sequence = "".join([sequence[cut_start_pos_list[i]:cut_end_pos_list[i]+1][::-1]
                                    if cut_start_pos_list[i] != cut_end_pos_list[i]
                                    else sequence[cut_start_pos_list[i]]
                                    for i in range(0, len(cut_start_pos_list))])
            output = ">rev_" + name + rev_sequence + "\n"
            return output
        # random shuffle the sequence
        case 3:
            name = input.split("\n")[0] + "\n"
            sequence = "".join(input.split("\n")[1:])
            shuff_sequence = "".join(random.sample(sequence,len(sequence)))
            output = ">rev_" + name + shuff_sequence + "\n"
            return output