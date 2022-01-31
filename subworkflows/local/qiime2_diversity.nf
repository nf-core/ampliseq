/*
 * Diversity indices with QIIME2
 */

include { QIIME2_TREE                 } from '../../modules/local/qiime2_tree'
include { QIIME2_ALPHARAREFACTION     } from '../../modules/local/qiime2_alphararefaction'
include { QIIME2_DIVERSITY_CORE       } from '../../modules/local/qiime2_diversity_core'
include { QIIME2_DIVERSITY_ALPHA      } from '../../modules/local/qiime2_diversity_alpha'
include { QIIME2_DIVERSITY_BETA       } from '../../modules/local/qiime2_diversity_beta'
include { QIIME2_DIVERSITY_ADONIS     } from '../../modules/local/qiime2_diversity_adonis'
include { QIIME2_DIVERSITY_BETAORD    } from '../../modules/local/qiime2_diversity_betaord'

workflow QIIME2_DIVERSITY {
    take:
    ch_metadata
    ch_asv
    ch_seq
    ch_stats //QIIME2_FILTERTAXA.out.tsv
    ch_metacolumn_pairwise //METADATA_PAIRWISE.out
    ch_metacolumn_all //METADATA_ALL.out
    skip_alpha_rarefaction
    skip_diversity_indices

    main:
    //Phylogenetic tree for beta & alpha diversities
    QIIME2_TREE ( ch_seq )

    //Alpha-rarefaction
    if (!skip_alpha_rarefaction) {
        QIIME2_ALPHARAREFACTION ( ch_metadata, ch_asv, QIIME2_TREE.out.qza, ch_stats )
    }

    //Calculate diversity indices
    if (!skip_diversity_indices) {

        QIIME2_DIVERSITY_CORE ( ch_metadata, ch_asv, QIIME2_TREE.out.qza, ch_stats )
        //Print warning if rarefaction depth is <10000
        QIIME2_DIVERSITY_CORE.out.depth.subscribe { if ( it.baseName.toString().startsWith("WARNING") ) log.warn it.baseName.toString().replace("WARNING ","QIIME2_DIVERSITY_CORE: ") }

        //alpha_diversity ( ch_metadata, DIVERSITY_CORE.out.qza, ch_metacolumn_all )
        ch_metadata
            .combine( QIIME2_DIVERSITY_CORE.out.vector.flatten() )
            .combine( ch_metacolumn_all )
            .set{ ch_to_diversity_alpha }
        QIIME2_DIVERSITY_ALPHA ( ch_to_diversity_alpha )

        //beta_diversity ( ch_metadata, DIVERSITY_CORE.out.qza, ch_metacolumn_pairwise )
        ch_metadata
            .combine( QIIME2_DIVERSITY_CORE.out.distance.flatten() )
            .combine( ch_metacolumn_pairwise )
            .set{ ch_to_diversity_beta }
        QIIME2_DIVERSITY_BETA ( ch_to_diversity_beta )

        //adonis ( ch_metadata, DIVERSITY_CORE.out.qza, ch_metacolumn_all )
        ch_metadata
            .combine( QIIME2_DIVERSITY_CORE.out.distance.flatten() )
            .combine( ch_metacolumn_pairwise )
            .set{ ch_to_diversity_beta }
        QIIME2_DIVERSITY_ADONIS ( ch_to_diversity_beta )

        //beta_diversity_ordination ( ch_metadata, DIVERSITY_CORE.out.qza )
        ch_metadata
            .combine( QIIME2_DIVERSITY_CORE.out.pcoa.flatten() )
            .set{ ch_to_diversity_betaord }
        QIIME2_DIVERSITY_BETAORD ( ch_to_diversity_betaord )
    }
}
