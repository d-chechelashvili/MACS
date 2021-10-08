//
//  CollectionViewCell.swift
//  Bidirectional Layout
//
//  Created by dato che on 12/16/20.
//

import UIKit

class CollectionViewAsCell: UITableViewCell {

@IBOutlet var collectionView: UICollectionView!

lazy var flowLayout: UICollectionViewFlowLayout = {
    let flowLayout = UICollectionViewFlowLayout()
    flowLayout.scrollDirection = .horizontal
    return flowLayout
}()


override func awakeFromNib() {
    super.awakeFromNib()
    configureCollectionView()
    // Initialization code
}

func configureCollectionView() {
    collectionView.delegate = self
    collectionView.dataSource = self
    
    collectionView.collectionViewLayout = flowLayout
    
    collectionView.register(
        UINib(nibName: "CollectionViewCellCell", bundle: nil),
        forCellWithReuseIdentifier: "CollectionViewCellCell"
    )
    
    collectionView.register(
        UINib(nibName: "CollectionViewCellCellRotated", bundle: nil),
        forCellWithReuseIdentifier: "CollectionViewCellCellRotated"
    )

}

override func setSelected(_ selected: Bool, animated: Bool) {
    super.setSelected(selected, animated: animated)

    // Configure the view for the selected state
}

}

extension CollectionViewAsCell: UICollectionViewDelegate {

}

extension CollectionViewAsCell: UICollectionViewDataSource {

func numberOfSections(in collectionView: UICollectionView) -> Int {
    return 1
}

func collectionView(_ collectionView: UICollectionView, numberOfItemsInSection section: Int) -> Int{
    return 6
}


func collectionView(_ collectionView: UICollectionView, cellForItemAt indexPath: IndexPath) -> UICollectionViewCell {
    if indexPath.row % 2 == 0{
        let cell = collectionView.dequeueReusableCell(withReuseIdentifier: "CollectionViewCellCell", for: indexPath)
        return cell
    }
    let cell = collectionView.dequeueReusableCell(withReuseIdentifier: "CollectionViewCellCellRotated", for: indexPath)
    return cell
}

}

extension CollectionViewAsCell: UICollectionViewDelegateFlowLayout {

func collectionView(
    _ collectionView: UICollectionView,
    layout collectionViewLayout: UICollectionViewLayout,
    sizeForItemAt indexPath: IndexPath
) -> CGSize {
    return CGSize(width: collectionView.frame.width/4, height: collectionView.frame.height)
}

}

