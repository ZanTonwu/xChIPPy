import pandas as pd
import pyliftover

def Converter(genome_build_in='mm9', genome_build_out='mm10'):
    return pyliftover.LiftOver(genome_build_in, genome_build_out)

def convertEnhancerDF(data, Converter, feature_name=None):
    cols=list(data.columns)
    if 'name' in data.columns:
        cols.remove('name')
    if len(cols) > 6:
        cols=cols[:6]
    cvrt_data = pd.DataFrame(columns=cols)
    cvrt_data.index.name='name'
    cvrt_idx = 1
    for idx in range(len(data)):
        cvrt_start = Converter.convert_coordinate(data.loc[idx]['chr'], data.loc[idx]['start'])
        cvrt_end = Converter.convert_coordinate(data.loc[idx]['chr'], data.loc[idx]['end'])
        if (cvrt_start is not None) and (cvrt_end is not None):
            if len(cvrt_start)==1 and len(cvrt_end)==1:                  # only one converted genomic region should exist
                if cvrt_start[0][0] == cvrt_end[0][0]:                  # feature must have both ends on the same chromosome
                    chrom = cvrt_start[0][0]
                    if cvrt_start[0][1] > cvrt_end[0][1]:
                        start = cvrt_end[0][1]
                        end = cvrt_start[0][1]
                    else:
                        start = cvrt_start[0][1]
                        end = cvrt_end[0][1]
                    if end > start:                                     # only take features that have a length greater than 0
                        if feature_name is not None:
                            index_name = feature_name + str(cvrt_idx)
                        else:
                            if 'name' in data.columns:
                                index_name=data.loc[idx]['name']
                            else:
                                index_name = cvrt_idx
                        row = [chrom, start, end]
                        if len(cols) > 3:
                            for col in cols[3:]:                        # transfer some information from columns
                                row.append(data.loc[idx][col])          # append new row with enhancer information, index is enhancer name
                        cvrt_data.loc[index_name]=row
                        cvrt_idx += 1
    return cvrt_data


enh_mm9 = pd.DataFrame.from_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/ES_Enhancers_mm9.csv", index_col=False)
enh_mm10 = convertEnhancerDF(enh_mm9, Converter('mm9', 'mm10'), feature_name='enh')
enh_mm10.to_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/ES_Enhancers_mm10.csv")

prL_mm9=pd.DataFrame.from_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/ES_PromoterLike_mm9.csv", index_col=False) 
prL_mm10 = convertEnhancerDF(prL_mm9, Converter('mm9', 'mm10'), feature_name='prL')
prL_mm10.to_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/ES_PromoterLike_mm10.csv")

