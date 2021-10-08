//
//  CollectionViewCell.swift
//  Contact Book
//
//  Created by dato che on 12/24/20.
//

import UIKit

class CollectionViewCell: UICollectionViewCell {
    
    @IBOutlet var initialsLabel: UILabel!
    @IBOutlet var numberLabel: UILabel!
    @IBOutlet var nameLabel: UILabel!
    @IBOutlet var wholeView: UIView!
    
    override func awakeFromNib() {
        super.awakeFromNib()
        // Initialization code
    }
    
    override var bounds: CGRect{
        didSet{
            layoutIfNeeded()
        }
    }

    override func layoutSubviews() {
        wholeView.layer.borderWidth = 1
        wholeView.layer.borderColor = UIColor.gray.cgColor
        wholeView.layer.cornerRadius = 8
        initialsLabel.layer.cornerRadius = initialsLabel.frame.width / 2
        initialsLabel.layer.masksToBounds = true
    }
    
}
