//
//  FirstSectionCell.swift
//  Bidirectional Layout
//
//  Created by dato che on 12/15/20.
//

import UIKit

class FirstSectionCell: UITableViewCell {
    
    @IBOutlet var view1: UIView!
    @IBOutlet var view2: UIView!

    override func awakeFromNib() {
        super.awakeFromNib()
        // Initialization code
        view1.layer.cornerRadius = 8
        view2.layer.cornerRadius = 8

    }

    override func setSelected(_ selected: Bool, animated: Bool) {
        super.setSelected(selected, animated: animated)

        // Configure the view for the selected state
    }
    
}
