import android.graphics.Rect
import android.view.View
import androidx.recyclerview.widget.RecyclerView


class SpacesItemDecoration(private val verticalSpace: Int, private val horizontalSpace: Int) :
    RecyclerView.ItemDecoration() {
    override fun getItemOffsets(
        outRect: Rect,
        view: View,
        parent: RecyclerView,
        state: RecyclerView.State
    ) {

        outRect.bottom = verticalSpace
        outRect.left = horizontalSpace / 2
        outRect.right = horizontalSpace / 2

        if (parent.getChildLayoutPosition(view) == 0 || parent.getChildLayoutPosition(view) == 1) {
            outRect.top = verticalSpace
        } else {
            outRect.top = 0
        }
    }
}