//
//  FullRedCell.swift
//  Bidirectional Layout
//
//  Created by dato che on 12/15/20.
//

import UIKit

class FullRedCell: UITableViewCell {

    @IBOutlet var redView: UIView!
    
    override func awakeFromNib() {
        super.awakeFromNib()
        redView.layer.cornerRadius = 6
        // Initialization code
    }

    override func setSelected(_ selected: Bool, animated: Bool) {
        super.setSelected(selected, animated: animated)

        // Configure the view for the selected state
    }
    
}
