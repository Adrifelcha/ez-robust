# MEGALEX Data

## Download Instructions

Download the following files from [lexique.org](http://www.lexique.org/?page_id=366):

1. **Megalex raw data.zip** (75 MB) - Contains trial-level data needed for computing RT variance/IQR
2. **Megalex-items-visual.tsv** - Item-level summary statistics

## Expected File Structure

After downloading and extracting, this folder should contain:

```
data/
├── README.md (this file)
├── Megalex-items-visual.tsv
└── megalex_visual_raw.csv (extracted from zip, may need renaming)
```

## Data Column Reference

### Item-level file (Megalex-items-visual.tsv)
Expected columns include:
- `word`: the stimulus word
- `rt`: mean response time
- `acc`: accuracy rate
- `freqlivres`: word frequency (books corpus)
- `nblettres`: number of letters

### Raw trial-level data
Expected columns:
- Participant identifier
- Item/word presented
- Response time (RT) in milliseconds
- Accuracy (1 = correct, 0 = incorrect)

## Notes

- The raw data is necessary for computing RT variance and IQR per condition
- Item-level summaries are useful for understanding the overall dataset structure
- Total dataset: ~28,466 words tested by ~100 participants each
