def check_quality(seqs, quality_threshold=0):
        new_seqs={}
        for name, (sequence, quality) in seqs.items():
            quality_num = []
            for symbol in quality:
                quality_num.append(ord(symbol)-33)
            mean_quality = sum(quality_num)/len(quality_num)
            if mean_quality >= quality_threshold:
                new_seqs[name]=(sequence, quality)
        return new_seqs



    def check_gc(new_sq, gc_bounds=(0, 100)):
        if not isinstance(gc_bounds, tuple):
            gc_bounds = (0, gc_bounds)
        seqs_after_gc={}
        for name, (sequence, quality) in new_sq.items():
            count = 0
            for nucleotide in sequence:
                if nucleotide in "CG":
                    count+=1
            percent = count/len(sequence)*100
            if gc_bounds[0] <= percent <= gc_bounds[1]:
                seqs_after_gc[name]=(sequence, quality)
        return seqs_after_gc





    def check_lenght(seqs_after_gc, length_bounds=(0, 2**32)):
        if not isinstance(length_bounds, tuple):
            length_bounds = (0, length_bounds)
        seqs_after_len={}
        for name, (sequence, quality) in seqs_after_gc.items():
            if length_bounds[0] <= len(sequence) <= length_bounds[1]:
                seqs_after_len[name]=(sequence, quality)
        return seqs_after_len


