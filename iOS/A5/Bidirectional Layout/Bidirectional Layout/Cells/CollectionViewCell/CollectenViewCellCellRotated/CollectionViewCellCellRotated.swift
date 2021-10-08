//
//  CollectionViewCellCellRotated.swift
//  Bidirectional Layout
//
//  Created by dato che on 12/15/20.
//

import UIKit

class CollectionViewCellCellRotated: UICollectionViewCell {

    
    @IBOutlet var smallerView: UIView!
    @IBOutlet var biggerView: UIView!
    
    override func awakeFromNib() {
        super.awakeFromNib()
        smallerView.layer.cornerRadius = 4
        biggerView.layer.cornerRadius = 4
        // Initialization code
    }

}
